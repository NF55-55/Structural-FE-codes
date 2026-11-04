#include "stiffnessMatrix.h"

namespace StiffMat{
    void computeKel(const NodesCoordEigen& coords, dNEigen& dN, invJEigen& invJ
        , BEigen& bEl, KelEigen& kEl)
    {
        kEl.setZero();
        Settings::precision dxdxi   {};
        Settings::precision dydxi   {};
        Settings::precision dxdeta  {};
        Settings::precision dydeta  {};
        Settings::precision detJ    {};
        Settings::precision invDetJ {};
        Settings::precision scalarValue {};
        for (const Settings::precision xi : gaussianPts){
            for (const Settings::precision eta : gaussianPts){
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
                scalarValue = Material::t*detJ; // Gauss weights omitted (=1)
                kEl.noalias() += scalarValue
                    *(bEl.transpose() * Material::C * bEl);
            }
        }
    }

    void computeKg(const Mesh& mesh, sparseMat& kG)
    {
        std::vector<Eigen::Triplet<Settings::precision>> tripletList;
        tripletList.reserve(mesh.nElements*Settings::nDofKel);
        StiffMat::KelEigen kEl;
        StiffMat::NodesCoordEigen nodesCoords;
        StiffMat::dNEigen dN;
        StiffMat::invJEigen invJ;
        StiffMat::BEigen bEl;
        Eigen::Matrix<int,Settings::dofNode,Settings::nNodesEl> 
            dofEliMat;
        int kElContiguousIndx {0};
        const float numEl1perc {mesh.nElements/20.0f};
        float numEl2Message {numEl1perc};
        for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
            nodesCoords = mesh.nodes(mesh.elements.row(e)
                ,Eigen::placeholders::all);
            computeKel(nodesCoords, dN, invJ, bEl, kEl);
            dofEliMat.row(0) = mesh.elements.row(e)*2;
            dofEliMat.row(1) = mesh.elements.row(e)*2 
                + Eigen::RowVector<int,Settings::nNodesEl>::Ones();
            const auto dofEli = dofEliMat.reshaped(); // view
            // as Eigen Matrices are column-major, stored contiguously 
            // column by column, it's more efficient to read them also 
            // in this way
            kElContiguousIndx = 0;
            for (const int j : dofEli){
                for (const int i : dofEli){
                    tripletList.emplace_back(i,j,kEl(kElContiguousIndx++));
                }
            }
            if (e > numEl2Message){
                std::cout << "Sparse Global Stiffness Matrix computation... " 
                    << numEl2Message/numEl1perc*5.0f << "%\n";
                numEl2Message += numEl1perc;
            }
        }
        kG.setFromTriplets(tripletList.begin(), tripletList.end());
    }
} 