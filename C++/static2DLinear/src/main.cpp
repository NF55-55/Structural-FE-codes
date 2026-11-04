#include <iostream>
#include <fstream>
#include "Eigen/Eigen"
#include "Eigen/PardisoSupport"
#include "Mesh.h"
#include "import.h"
#include "settings.h"
#include "TimerHighRes.h"
#include "estimateMemoryUsage.h"
#include "material.h"
#include "Bc.h"
#include "rearrangeDof.h"
#include "stiffnessMatrix.h"
#include "export.h"

int main(){
    //Timer timer {};
    
    Mesh mesh {};
    Bc bc {}; // boundary conditions
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    ndofVectf disp;
    ndofVectf force;
    
    // To handle exceptions in any function/method
    try{
        //-----Importing mesh and boundary conditions
        importMesh(mesh, bc);
        std::cout << "Problem size: " << mesh.nDof << " d.o.f.\n";
        double estimatedMemoryUsageMB {estimateMemUsageMB(mesh.nDof
            , mesh.nElements)};
        std::cout << "Esitmated Memory usage (before solving phase): " 
            << estimatedMemoryUsageMB/Settings::multiplierMemory4Solving 
            << " MB\nEstimated MAX Memory usage: " << estimatedMemoryUsageMB 
            << " MB\n\n";
        if (estimatedMemoryUsageMB > Settings::memoryLimitMB)
            throw "Estimated memory usage > memory usage limit set in settings";
        disp = ndofVectf::Constant(mesh.nDof,Settings::NaN);
        force = ndofVectf::Zero(mesh.nDof);

        bc.setBc(disp, force);
        
        // Partitioning method used for solving linear system of equations.
        // Re-arranging nodal disp. vector (unknown at the top and known vals 
        // at the bottom) and nodal force vector (opposite arrangement) + saving
        // the mapping vector (equal for both vectors)
        ndofVecti mapDof;
        mapDof = ndofVecti::LinSpaced(mesh.nDof, 0, mesh.nDof-1);
        rearrange(disp, force, mapDof);

        //-----Sparse Global stiffness matrix computation + solving
        Eigen::SparseMatrix<Settings::precision> kG (mesh.nDof, mesh.nDof);
        StiffMat::computeKg(mesh, kG);
        
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> 
            mapDofMat {mapDof};
        kG = kG.twistedBy(mapDofMat.transpose()); // Re-arranging
        std::cout << "\nGlobal Stiffness Matrix [K] computed "
            << "successfully!\n\n";

        //-----Solving linear system of eqs. (partitioning method)
        int lastNanIndx {mesh.nDof-1};
        for (; lastNanIndx>0; --lastNanIndx){
            if (std::isnan(disp(lastNanIndx))) break;
        }
        int nUnknownDisp {lastNanIndx+1};
        int nKnownDisp {mesh.nDof - nUnknownDisp};
        std::cout << "Solving the system of linear algebraic equations [K]"
            << "{u}={F}...\n\n";
        
        Settings::solverType solver;
        solver.compute(kG.block(0,0,lastNanIndx+1,lastNanIndx+1));
        if(solver.info()!=Eigen::Success) {
            throw "Linear solver failed decomposition(1).";
        }
        Eigen::Vector<Settings::precision,Eigen::Dynamic> rhsDispSol = 
            force.head(nUnknownDisp);
        rhsDispSol.noalias() -= kG.topRightCorner(nUnknownDisp,nKnownDisp)
            *disp.tail(nKnownDisp);
        disp.head(nUnknownDisp) = solver.solve(rhsDispSol);
        if(solver.info()!=Eigen::Success) {
            throw "Linear solver failed solving(1).";
        }
        std::cout << "Unknown nodal displacements computed successfully!"
            << "\n\n";
        disp = mapDofMat * disp;

        //-----Exporting the results
        if (Settings::exportBinaryVTK) Export::exportResultsBinary(mesh, disp);
        else Export::exportResultsASCII(mesh, disp);
    }
    catch(const char* exceptionMsg){
        std::cerr << "Error: " << exceptionMsg << '\n';
    }
    catch(...){
        std::cerr << "Unrecognized error occured!\n";
    }

    return 0;
}
