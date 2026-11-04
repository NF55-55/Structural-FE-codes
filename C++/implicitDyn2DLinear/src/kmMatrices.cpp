#include "kmMatrices.h"

namespace kmMat{
    void computeElKM(const NodesCoordEigen& coords, NEigen& N, dNEigen& dN
        , invJEigen& invJ, BEigen& bEl, KMelEigen& kEl, KMelEigen& mEl)
    {
        kEl.setZero();
        mEl.setZero();
        Settings::precision dxdxi   {};
        Settings::precision dydxi   {};
        Settings::precision dxdeta  {};
        Settings::precision dydeta  {};
        Settings::precision detJ    {};
        Settings::precision invDetJ {};
        Settings::precision scalarValueK {};
        Settings::precision scalarValueM {};
        for (const Settings::precision xi : gaussianPts){
            for (const Settings::precision eta : gaussianPts){
                N << (1.0-xi)*(1.0-eta),0.0,(1.0+xi)*(1.0-eta),0.0
                    ,(1.0+xi)*(1.0+eta),0.0,(1.0-xi)*(1+eta),0.0
                    ,0.0,(1.0-xi)*(1.0-eta),0.0,(1.0+xi)*(1.0-eta)
                    ,0.0,(1.0+xi)*(1.0+eta),0.0,(1.0-xi)*(1+eta);
                N *= 0.25;
                dN <<(-1.0+eta),0.0,(1.0-eta),0.0,(1.0+eta),0.0,(-1.0-eta),0.0,
                    (-1.0+xi),0.0,(-1.0-xi),0.0,(1.0+xi),0.0,(1.0-xi),0.0,
                    0.0,(-1.0+eta),0.0,(1.0-eta),0.0,(1.0+eta),0.0,(-1.0-eta),
                    0.0,(-1.0+xi),0.0,(-1.0-xi),0.0,(1.0+xi),0.0,(1.0-xi);
                dN *= 0.25;
                dxdxi = 0.25*(-coords(0,0)+coords(1,0)+coords(2,0)-coords(3,0)
                    +eta*(coords(0,0)-coords(1,0)+coords(2,0)-coords(3,0)));
                dydxi = 0.25*(-coords(0,1)+coords(1,1)+coords(2,1)-coords(3,1)
                    +eta*(coords(0,1)-coords(1,1)+coords(2,1)-coords(3,1)));
                dxdeta = 0.25*(-coords(0,0)-coords(1,0)+coords(2,0)+coords(3,0)
                    +xi*(coords(0,0)-coords(1,0)+coords(2,0)-coords(3,0)));
                dydeta = 0.25*(-coords(0,1)-coords(1,1)+coords(2,1)+coords(3,1)
                    +xi*(coords(0,1)-coords(1,1)+coords(2,1)-coords(3,1)));
                detJ = dxdxi*dydeta - dydxi*dxdeta;
                invDetJ = 1.0/detJ;
                invJ << dydeta, -dydxi, 0, 0,
                        -dxdeta, dxdxi, 0, 0,
                        0, 0, dydeta, -dydxi,
                        0, 0, -dxdeta, dxdxi;
                invJ *= invDetJ;
                bEl.noalias() = arrangeMat * invJ * dN;
                scalarValueK = Material::t*detJ; // Gauss weights omitted (=1)
                scalarValueM = Material::t*detJ*Material::rho; // same here
                kEl.noalias() += scalarValueK
                    *(bEl.transpose() * Material::C * bEl);
                mEl.noalias() += scalarValueM * (N.transpose() * N);
            }
        }
    }

    void computekGmG(const Mesh& mesh, sparseMat& mG, sparseMat& kG)
    {
        std::vector<Eigen::Triplet<Settings::precision>> tripletListK;
        std::vector<Eigen::Triplet<Settings::precision>> tripletListM;
        tripletListK.reserve(mesh.nElements*Settings::nDofKel);
        tripletListM.reserve(mesh.nElements*Settings::nDofKel);
        kmMat::KMelEigen kEl;
        kmMat::KMelEigen mEl;
        kmMat::NodesCoordEigen nodesCoords;
        kmMat::NEigen N;
        kmMat::dNEigen dN;
        kmMat::invJEigen invJ;
        kmMat::BEigen bEl;
        Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
            dofEliMat;
        int kmElContiguousIndx {};
        const float numEl1perc {mesh.nElements/20.0f};
        float numEl2Message {numEl1perc};
        for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
            nodesCoords = mesh.nodes(mesh.elements.row(e)
                ,Eigen::placeholders::all);
            kmMat::computeElKM(nodesCoords,N,dN,invJ,bEl,kEl,mEl);
            dofEliMat.row(0) = mesh.elements.row(e)*2;
            dofEliMat.row(1) = mesh.elements.row(e)*2 
                + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
            const auto dofEli = dofEliMat.reshaped(); // view
            // as Eigen Matrices are column-major, stored contiguously 
            // column by column, it's more efficient to read them also 
            // in this way
            kmElContiguousIndx = 0;
            for (const int j : dofEli){
                for (const int i : dofEli){
                    tripletListM.emplace_back(i,j,mEl(kmElContiguousIndx));
                    tripletListK.emplace_back(i,j,kEl(kmElContiguousIndx));
                    ++kmElContiguousIndx;
                }
            }
            if (e > numEl2Message){
                std::cout << "Sparse global mass & stiffness matrix " 
                    << "computation... " << numEl2Message/numEl1perc
                    *5.0f << "%\n";
                numEl2Message += numEl1perc;
            }
        }
        kG.setFromTriplets(tripletListK.begin(), tripletListK.end());
        mG.setFromTriplets(tripletListM.begin(), tripletListM.end());
        std::cout << "\nGlobal Mass [M] and Stiffness [K] matrices computed"
            << " successfully!\n\n";
    }
} 