#ifndef STIFFNESS_MATRIX_H
#define STIFFNESS_MATRIX_H

#include <array>
#include <iostream>
#include "Eigen/Core"
#include "Mesh.h"
#include "settings.h"
#include "material.h"

namespace StiffMat{
    using SizT = std::size_t;
    using sparseMat = Eigen::SparseMatrix<Settings::precision>;
    using KelEigen = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl>;
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
    void computeKg(const Mesh& mesh, sparseMat& kG);
}


#endif