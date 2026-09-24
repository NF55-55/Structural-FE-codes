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
#include "elStiffnessMatrix.h"
#include "export.h"

// to set number of threads for intel mkl pardiso sparse linear solver
//extern "C" void mkl_set_num_threads(int nth);

int main(){
    //mkl_set_num_threads(DEFAULT_MKL_THREADS);
    
    //Timer timer {};
    //timer.start();
    Mesh mesh {};
    Bc bc {}; // boundary conditions
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using SizT = std::size_t;
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

        // Re-arranging nodal disp. vector (unknown at the top and known vals 
        // at the bottom) and nodal force vector (opposite arrangement) + saving
        // the mapping vector (equal for both vectors)
        ndofVecti mapDof;
        mapDof = ndofVecti::LinSpaced(mesh.nDof, 0, mesh.nDof-1);
        rearrange(disp, force, mapDof);

        //-----Sparse Global stiffness matrix computation + solving
        {
            Eigen::SparseMatrix<Settings::precision> kG (mesh.nDof, mesh.nDof);
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
                const float numEl1perc {mesh.nElements/50.0f};
                float numEl2Message {numEl1perc};
                for (SizT e=0uz; e<static_cast<SizT>(mesh.nElements); ++e){
                    nodesCoords = mesh.nodes(mesh.elements.row(e)
                        ,Eigen::placeholders::all);
                    StiffMat::computeKel(nodesCoords, dN, invJ, bEl, kEl);
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
                            tripletList.emplace_back(i,j
                                ,kEl(kElContiguousIndx++));
                        }
                    }
                    if (e > numEl2Message){
                        std::cout << "Sparse Global Stiffness Matrix " 
                            << "computation... " << numEl2Message/numEl1perc
                            *2.0f << "\%\n";
                        numEl2Message += numEl1perc;
                    }
                }
                kG.setFromTriplets(tripletList.begin(), tripletList.end());
            } // Triplet list and all other helper matrices deallocated
            Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> 
                mapDofMat {mapDof};
            kG = kG.twistedBy(mapDofMat.transpose()); // Re-arranging
            std::cout << "\nGlobal Stiffness Matrix [K] computed "
                << "successfully!\n\n";

            //-----Solving linear system of eqs. (splitting method)
            int lastNanIndx {mesh.nDof-1};
            for (; lastNanIndx>0; --lastNanIndx){
                if (std::isnan(disp(lastNanIndx))) break;
            }
            int nUnknownDisp {lastNanIndx+1};
            int nKnownDisp {mesh.nDof - nUnknownDisp};
            std::cout << "Solving the system of linear algebraic equations [K]"
                << "{u}={F}...\n\n";
            // Suggested for square positive definite system:
            // For small problems (up to 200k dof) -> SimplicialLLT
            // For medium problems (up to 1M dof)  -> PardisoLLT (linking)
            // For large problems (over 1-2M dof)  -> ConjugateGradient
            Eigen::PardisoLLT<Eigen::SparseMatrix<Settings::precision>> 
                solver;
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
        } // kG and some other matrices deallocated

        //-----Exporting the results
        if (Settings::exportBinaryVTK) Export::exportResultsBinary(mesh, disp);
        else Export::exportResultsASCII(mesh, disp);
        
        //timer.stop();
        //std::cout << "Elapsed: " << std::setprecision(10) 
        //    << timer.getDurationMicro() << " micro sec.\n";
    }
    catch(const char* exceptionMsg){
        std::cerr << "Error: " << exceptionMsg << '\n';
    }
    catch(...){
        std::cerr << "Unrecognized error occured!\n";
    }

    return 0;
}
