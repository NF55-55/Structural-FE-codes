#ifndef EL_STIFFNESS_MATRIX_H
#define EL_STIFFNESS_MATRIX_H

#include <array>
#include "Eigen/Core"
#include "settings.h"
#include "material.h"

namespace StiffMat{
    using KelEigen = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl,Eigen::RowMajor>;
    using NodesCoordEigen = Eigen::Matrix<Settings::precision
        ,Settings::nNodesEl,Settings::dofNode>;
    using dNEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nDofEl>;
    using invJEigen = Eigen::Matrix<Settings::precision,Settings::nNodesEl
        ,Settings::nNodesEl>;
    using BEigen = Eigen::Matrix<Settings::precision,3uz,Settings::nDofEl>;

    constexpr const std::array<Settings::precision,2uz> gaussianPts 
        {{-1.0/Settings::sqrt3, 1.0/Settings::sqrt3}};

    constexpr const Eigen::Matrix<Settings::precision,3uz,4uz> arrangeMat 
    {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0, 1.0},
        {0.0, 1.0, 1.0, 0.0}
    };

    void computeKel(const NodesCoordEigen& coords, dNEigen& dN, invJEigen& invJ
        , BEigen& bEl, KelEigen& kEl);
}


#endif