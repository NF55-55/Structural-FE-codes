#include "massMatrix.h"

namespace mMat{
    void computeElM(const NodesCoordEigen& coords, NEigen& N, MelEigen& mEl)
    {
        mEl.setZero();
        Settings::precision dxdxi   {};
        Settings::precision dydxi   {};
        Settings::precision dxdeta  {};
        Settings::precision dydeta  {};
        Settings::precision detJ    {};
        Settings::precision scalarValueM {};
        for (const Settings::precision xi : gaussianPts){
            for (const Settings::precision eta : gaussianPts){
                N << (1.0-xi)*(1.0-eta),0.0,(1.0+xi)*(1.0-eta),0.0
                    ,(1.0+xi)*(1.0+eta),0.0,(1.0-xi)*(1+eta),0.0
                    ,0.0,(1.0-xi)*(1.0-eta),0.0,(1.0+xi)*(1.0-eta)
                    ,0.0,(1.0+xi)*(1.0+eta),0.0,(1.0-xi)*(1+eta);
                N *= 0.25;
                dxdxi = 0.25*(-coords(0,0)+coords(1,0)+coords(2,0)-coords(3,0)
                    +eta*(coords(0,0)-coords(1,0)+coords(2,0)-coords(3,0)));
                dydxi = 0.25*(-coords(0,1)+coords(1,1)+coords(2,1)-coords(3,1)
                    +eta*(coords(0,1)-coords(1,1)+coords(2,1)-coords(3,1)));
                dxdeta = 0.25*(-coords(0,0)-coords(1,0)+coords(2,0)+coords(3,0)
                    +xi*(coords(0,0)-coords(1,0)+coords(2,0)-coords(3,0)));
                dydeta = 0.25*(-coords(0,1)-coords(1,1)+coords(2,1)+coords(3,1)
                    +xi*(coords(0,1)-coords(1,1)+coords(2,1)-coords(3,1)));
                detJ = dxdxi*dydeta - dydxi*dxdeta;
                scalarValueM = Material::t*detJ*Material::rho; // gauss w = 1
                mEl.noalias() += scalarValueM * (N.transpose() * N);
            }
        }
    }
    void computeGlobalM(SparseMat& mG, const Mesh& mesh)
    {
        std::vector<Eigen::Triplet<Settings::precision>> tripletListM;
        tripletListM.reserve(mesh.nElements*Settings::nDofKel);
        MelEigen mEl;
        NodesCoordEigen nodesCoords;
        NEigen N;
        Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
            dofEliMat;
        int mElContiguousIndx {};
        for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
            nodesCoords = mesh.nodes(mesh.elements.row(e)
                ,Eigen::placeholders::all);
            computeElM(nodesCoords,N,mEl);
            dofEliMat.row(0) = mesh.elements.row(e)*2;
            dofEliMat.row(1) = mesh.elements.row(e)*2 
                + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
            const auto dofEli = dofEliMat.reshaped(); // view
            // as Eigen Matrices are column-major, stored contiguously 
            // column by column, it's more efficient to read them also 
            // in this way
            mElContiguousIndx = 0;
            for (const int j : dofEli){
                for (const int i : dofEli){
                    tripletListM.emplace_back(i,j
                        ,mEl(mElContiguousIndx++));
                }
            }
        }
        mG.setFromTriplets(tripletListM.begin(), tripletListM.end());
        std::cout << "\nGlobal Mass matrix [M] computed successfully!\n\n";
    }
} 