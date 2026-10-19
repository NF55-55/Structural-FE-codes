#include "elStiffnessMatrix.h"

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
} 