#ifndef SETTINGS_H
#define SETTINGS_H

#include <limits>
#include <filesystem>
#include "Eigen/Eigen"
#include "Eigen/PardisoSupport"

namespace Settings{
    // some aliases
    using precision = double; // for computations (float or double), leave
    // double when using Intel MKL Pardiso solvers
    using fPath = std::filesystem::path;
    // Suggested for square positive definite system:
    // For problems up to 1M dof -> PardisoLLT (direct)
    // For large problems (over 1-2M dof) -> ConjugateGradient (iterative)
    using solverType = Eigen::PardisoLLT<Eigen::SparseMatrix<precision>>;

    // export options and more
    constexpr const bool exportBinaryVTK {true}; // binary is much faster
    constexpr const bool export4dynamics {false};
    constexpr const bool fullPrecisionASCIIWriting {false}; // suggested false
    constexpr const double memoryLimitMB {3000.0}; // avoiding dynamic 
    // allocation errors due to lack of heap memory availability
    inline fPath expectedProjRoot {"static2DLinear"};
    inline fPath disp4DynPath {"./static2dynamic/nodalDispVect"};

    // Useful constants:
    constexpr const double multiplierMemory4Solving {2.5}; // the multiplier of 
    // memory estimated before solving, to account for the one of solving phase
    constexpr const precision NaN {std::numeric_limits<precision>::quiet_NaN()};
    constexpr const precision sqrt3 {1.73205080756887729352};
    constexpr const std::size_t nNodesEl {4uz};
    constexpr const std::size_t dofNode {2uz};
    constexpr const std::size_t nDofEl {nNodesEl * dofNode};
    constexpr const std::size_t nDofKel {nDofEl * nDofEl};
}

#endif