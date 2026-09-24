#ifndef SETTINGS_H
#define SETTINGS_H

#include <iostream>
#include <sstream>
#include <omp.h>
#include <limits>
#include <Eigen/Core>

namespace Settings{
    using precision = double; // for computations (float or double), leave
    // double when using Intel MKL Pardiso solvers
    constexpr const bool exportBinaryVTK{true}; // binary is much faster
    constexpr const bool fullPrecisionASCIIWriting {false}; // suggested false
    constexpr const double memoryLimitMB {3000.0}; // avoiding dynamic 
    // allocation errors due to lack of heap memory availability

    // Useful constants:
    constexpr const double multiplierMemory4Solving {2.5}; // the multiplier of 
    // memory estimated before solving, to account for the one of solving phase
    constexpr const precision NaN {std::numeric_limits<precision>::quiet_NaN()};
    constexpr const precision sqrt3 {1.73205080756887729352};
    constexpr const std::size_t nNodesEl {4uz};
    constexpr const std::size_t dofNode {2uz};
    constexpr const std::size_t nDofEl {nNodesEl * dofNode};
    constexpr const std::size_t nDofKel {nDofEl * nDofEl};

    void setNumThreads(int argc, char* argv[]);
}

#endif