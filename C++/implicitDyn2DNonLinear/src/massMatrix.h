#ifndef MASS_MATRIX_H
#define MASS_MATRIX_H

#include <array>
#include <vector>
#include <iostream>
#include "Eigen/Core"
#include "settings.h"
#include "material.h"
#include "Mesh.h"

namespace mMat{
    using SizT = std::size_t;
    using MelEigen = Eigen::Matrix<Settings::precision,Settings::nDofEl
        ,Settings::nDofEl>;
    using NodesCoordEigen = Eigen::Matrix<Settings::precision
        ,Settings::nNodesEl,Settings::dofNode>;
    using NEigen = Eigen::Matrix<Settings::precision,Settings::dofNode
        ,Settings::nDofEl>;
    using SparseMat = Eigen::SparseMatrix<Settings::precision>;

    constexpr const std::array<Settings::precision,2uz> gaussianPts 
        {{-1.0/Settings::sqrt3, 1.0/Settings::sqrt3}};

    void computeElM(const NodesCoordEigen& coords, NEigen& N, MelEigen& mEl);
    void computeGlobalM(SparseMat& mG, const Mesh& mesh);
}

#endif