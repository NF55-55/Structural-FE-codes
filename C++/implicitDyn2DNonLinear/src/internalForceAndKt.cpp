#include "internalForceAndKt.h"

namespace FintKtInt{
    void KtIntCache::checkAndResetValues(const sparseMat& ktInt
        , const Mesh& mesh)
    {
        bool isktIntCompatible = ktInt.isCompressed() && (ktInt.rows()==m_rows)
            && (ktInt.cols()==m_cols) && (ktInt.nonZeros()==m_nnz) 
            && (ktInt.innerIndexPtr()==m_innerIndices) 
            && (ktInt.outerIndexPtr()==m_outerIndices)
            && (ktInt.valuePtr()==m_values);
        if (isktIntCompatible){
            for (sparseMat::StorageIndex i=0; i<m_nnz; ++i){
                *(m_values + i) = 0.0;
            }
            if (isCached() && !m_isLoopInnerIndxCached 
                && m_cachedLoopInnerIndx.empty())
            {
                m_cachedLoopInnerIndx.resize(static_cast<SizT>(mesh.nElements)
                    *Settings::nDofKel);
                m_nElements = mesh.nElements;
            }
        }else throw "ktInt sparse matrix has changed, cannot use cache!";
    }

    void KtIntCache::addContribution2El(
        const sparseMat::StorageIndex row, const sparseMat::StorageIndex col
        , const SizT elIndx, const sparseMat::StorageIndex localContigIndx
        , const Settings::precision val)
    {
        using StoIndx = sparseMat::StorageIndex;
        SizT localIndxSizt {static_cast<SizT>(localContigIndx)};
        if (m_isLoopInnerIndxCached){
            StoIndx foundCachedInnerIndx {m_cachedLoopInnerIndx[elIndx
                *Settings::nDofKel + localContigIndx]};
        #pragma omp atomic update
            *(m_values + foundCachedInnerIndx) += val;
        }else{
            StoIndx nnzUpToSelectedCol {*(m_outerIndices + col)};
            StoIndx nnzUpToNextCol {*(m_outerIndices + col + 1)};

            StoIndx* first {m_innerIndices + nnzUpToSelectedCol};
            StoIndx* last  {m_innerIndices + nnzUpToNextCol}; // pointing to 
            // first nnz of next column
            
            StoIndx* inner {std::lower_bound(first,last,row)}; // binary search
            if (inner != last){
                StoIndx foundIndx {static_cast<StoIndx>(inner-m_innerIndices)};
                m_cachedLoopInnerIndx[elIndx*Settings::nDofKel+localIndxSizt] 
                    = foundIndx;
            #pragma omp atomic update
                *(m_values + foundIndx) += val;
            } else throw "Element not found in cache!";
        }
    }

    void compute(const Mesh& mesh, const nDofVectf& dispNew, nDofVectf& fInt
        , sparseMat& ktInt, KtIntCache& ktIntCache)
    {
        bool isKtIntCached {ktIntCache.isCached()};
        static const SizT elPerThread {static_cast<SizT>(
            mesh.nElements/Settings::nThreads)+1uz};
        static const DNprecomputed s_dN {};
        std::vector<nDofVectf> threadFintVect (
            static_cast<SizT>(Settings::nThreads));
        std::vector<std::vector<Eigen::Triplet<Settings::precision>>> 
            threadTripletKtIntVect (Settings::nThreads); // unused with cache

        for (SizT threadID=0uz; threadID<Settings::nThreads; ++threadID){
            threadFintVect[threadID] = nDofVectf::Zero(fInt.size());
            if (!isKtIntCached){
                threadTripletKtIntVect[threadID].reserve(elPerThread
                    *Settings::nDofKel);
            }
        }
        fInt.setZero();
        if (isKtIntCached) ktIntCache.checkAndResetValues(ktInt, mesh);
        // estimated memory usage of local variables except tripletlist 
        // (dynamically allocated), for each thread: 1052 B
    #pragma omp parallel
        {
        SizT threadID {static_cast<SizT>(omp_get_thread_num())};
        nDofVectf& threadFint = threadFintVect[threadID];
        auto& threadTripletList = threadTripletKtIntVect[threadID];
        Settings::precision dXdxi   {};
        Settings::precision dYdxi   {};
        Settings::precision dXdeta  {};
        Settings::precision dYdeta  {};
        Settings::precision detJ    {};
        Settings::precision invDetJ {};
        invJEigen invJ;
        B0Eigen b0mat;
        voigtTensorFull fVect;
        tensor fTensor; // [F]
        tensor eTensor; // [E]
        voigtTensor eVect;
        tensor pk2Tensor; // 2nd Piola-Kirchhoff stress tensor
        voigtTensor pk2Vect;
        tensor pTensor; // Nominal stress tensor
        voigtTensorFull pVect;
        nDofElVectf b0p;
        Settings::precision scalar4GaussInt {};
        matMultipliedByVel eDotStar;
        matMultipliedByVel pk2DotStar;
        matMultipliedByVelFull pDotStar;
        matMultipliedByVelFull fDotStar;
        elMat ktMat, ktGeo;
        mMat::NodesCoordEigen coordsL;
        nNodesRowVecti nodesID;
        nNodesRowVecti onesRowVect = nNodesRowVecti::Ones();
        nDofElVectf coordsE;
        Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
            dofEliMat;
        SizT dNindx {};
        int elMatContiguousIndx {};

    #pragma omp for schedule (static)
        for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
            ktMat.setZero();
            ktGeo.setZero();
            // Lagrangian coordinates (upper case X and Y in the notes)
            nodesID = mesh.elements.row(e);
            coordsL = mesh.nodes(nodesID,Eigen::placeholders::all);
            dofEliMat.row(0).noalias() = nodesID*Settings::dofNode;
            dofEliMat.row(1).noalias() = nodesID*Settings::dofNode+onesRowVect;
            const auto dofEli = dofEliMat.reshaped(); // view
            // Eulerian coordinates (lower case x and y in the notes)
            coordsE = coordsL.transpose().reshaped() + dispNew(dofEli);
            dNindx = 0uz;
            
            for (const auto xi : mMat::gaussianPts){
                for (const auto eta : mMat::gaussianPts){
                    dXdxi = 0.25*(-coordsL(0,0)+coordsL(1,0)+coordsL(2,0)
                        -coordsL(3,0)+eta*(coordsL(0,0)-coordsL(1,0)
                        +coordsL(2,0)-coordsL(3,0)));
                    dYdxi = 0.25*(-coordsL(0,1)+coordsL(1,1)+coordsL(2,1)
                        -coordsL(3,1)+eta*(coordsL(0,1)-coordsL(1,1)
                        +coordsL(2,1)-coordsL(3,1)));
                    dXdeta = 0.25*(-coordsL(0,0)-coordsL(1,0)+coordsL(2,0)
                        +coordsL(3,0)+xi*(coordsL(0,0)-coordsL(1,0)+coordsL(2,0)
                        -coordsL(3,0)));
                    dYdeta = 0.25*(-coordsL(0,1)-coordsL(1,1)+coordsL(2,1)
                        +coordsL(3,1)+xi*(coordsL(0,1)-coordsL(1,1)+coordsL(2,1)
                        -coordsL(3,1)));
                    detJ = dXdxi*dYdeta - dYdxi*dXdeta;
                    invDetJ = 1.0/detJ;
                    invJ << dYdeta, -dYdxi, 0, 0,
                            -dXdeta, dXdxi, 0, 0,
                            0, 0, dYdeta, -dYdxi,
                            0, 0, -dXdeta, dXdxi;
                    invJ *= invDetJ;
                    b0mat.noalias() = invJ * s_dN.getdNi(dNindx++);
                    fVect.noalias() = b0mat * coordsE;
                    fTensor << fVect(0),fVect(1),fVect(2),fVect(3);
                    eTensor.noalias() = 0.5 * fTensor.transpose() * fTensor;
                    eTensor.noalias() -= 0.5 
                        * tensor::Identity(Settings::dofNode,Settings::dofNode);
                    eVect << eTensor(0,0), eTensor(1,1), 2.0*eTensor(0,1);
                    pk2Vect.noalias() = Material::C * eVect;
                    pk2Tensor << pk2Vect(0), pk2Vect(2), pk2Vect(2), pk2Vect(1);
                    pTensor.noalias() = fTensor * pk2Tensor;
                    pVect <<pTensor(0,0),pTensor(0,1),pTensor(1,0),pTensor(1,1);
                    scalar4GaussInt = Material::t * detJ; // gauss weights = 1
                    b0p.noalias() = b0mat.transpose() * pVect;
                    // computing F_int for this element
                    for (int i=0; i<Settings::nDofEl; ++i){
                        threadFint(dofEli(i)) += scalar4GaussInt * b0p(i);
                    }
                    // computing material tangent stiffness matrix
                    eDotStar.row(0) = b0mat.row(0)*fTensor(0,0) + b0mat.row(2)
                        *fTensor(1,0);
                    eDotStar.row(1) = b0mat.row(1)*fTensor(0,1) + b0mat.row(3)
                        *fTensor(1,1);
                    eDotStar.row(2) = b0mat.row(0)*fTensor(0,1) + b0mat.row(2)
                        *fTensor(1,1) + b0mat.row(1)*fTensor(0,0) +b0mat.row(3)
                        *fTensor(1,0);
                    pk2DotStar.noalias() = Material::C * eDotStar;
                    pDotStar.row(0) = pk2DotStar.row(0)*fTensor(0,0) 
                        + pk2DotStar.row(2)*fTensor(0,1);
                    pDotStar.row(1) = pk2DotStar.row(2)*fTensor(0,0)
                        + pk2DotStar.row(1)*fTensor(0,1);
                    pDotStar.row(2) = pk2DotStar.row(0)*fTensor(1,0)
                        + pk2DotStar.row(2)*fTensor(1,1);
                    pDotStar.row(3) = pk2DotStar.row(2)*fTensor(1,0)
                        + pk2DotStar.row(1)*fTensor(1,1);
                    ktMat.noalias() += scalar4GaussInt * b0mat.transpose() 
                        * pDotStar;
                    // computing geometric tangent stiffness matrix
                    fDotStar.row(0) = b0mat.row(0)*pk2Tensor(0,0) 
                        + b0mat.row(1)*pk2Tensor(1,0);
                    fDotStar.row(1) = b0mat.row(0)*pk2Tensor(0,1)
                        + b0mat.row(1)*pk2Tensor(1,1);
                    fDotStar.row(2) = b0mat.row(2)*pk2Tensor(0,0)
                        + b0mat.row(3)*pk2Tensor(1,0);
                    fDotStar.row(3) = b0mat.row(2)*pk2Tensor(0,1)
                        + b0mat.row(3)*pk2Tensor(1,1);
                    ktGeo.noalias() += scalar4GaussInt * b0mat.transpose()
                        * fDotStar;
                }
            }
            
            elMatContiguousIndx = 0;
            for (const auto j : dofEli){
                for (const auto i : dofEli){
                    if (isKtIntCached){
                        ktIntCache.addContribution2El(i,j
                        ,e,elMatContiguousIndx
                        ,ktMat(elMatContiguousIndx)+ktGeo(elMatContiguousIndx));
                    }else{
                        threadTripletList.emplace_back(i,j
                        ,ktMat(elMatContiguousIndx)+ktGeo(elMatContiguousIndx));
                    }
                    ++elMatContiguousIndx;
                }
            }
        }
    #pragma omp for schedule (static) 
        for (int dof=0; dof<fInt.size(); dof++){
            // avoiding race condition by distributing dof to threads and not 
            // fInt local(thread) vectors
            Settings::precision forceSum {};
            for (SizT thread=0uz; thread<Settings::nThreads; ++thread){
                forceSum += threadFintVect[thread](dof);
            }
            fInt(dof) = forceSum;
        }
        }
        
        if (ktIntCache.isCached() && !ktIntCache.isInnerIndxCached()){
            ktIntCache.setInnerIndxToCached();
        }

        if (!isKtIntCached){
            std::vector<Eigen::Triplet<Settings::precision>> tripletListKtInt;
            tripletListKtInt.reserve(mesh.nElements*Settings::nDofKel);
            for (const auto& threadTriplets : threadTripletKtIntVect){
                tripletListKtInt.insert(tripletListKtInt.end()
                    , threadTriplets.begin(), threadTriplets.end());
            }
            ktInt.setFromTriplets(tripletListKtInt.begin()
                , tripletListKtInt.end());
            if (!ktInt.isCompressed()) ktInt.makeCompressed();

            ktIntCache.setupCache(ktInt);
        }
    }
}