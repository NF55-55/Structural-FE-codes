#include <iostream>
#include <string_view>
#include <algorithm>
#include "Eigen/Eigen"
#include "Eigen/PardisoSupport"
#include "Mesh.h"
#include "import.h"
#include "settings.h"
#include "TimerHighRes.h"
#include "material.h"
#include "Bc.h"
#include "massMatrix.h"
#include "export.h"
#include "Ic.h"
#include "clearResults.h"
#include "internalForceAndKt.h"

int main(int argc, char* argv[]){
    //Timer timer {};

    Mesh mesh {};
    Bc bc {}; // boundary conditions class
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using SizT = std::size_t;
    
    // To handle exceptions in any function/method
    try{
        // OpenMP multi-thread settings
        Settings::setNumThreadsFromUserInput(argc, argv);
        //-----Importing mesh and boundary conditions
        importMesh(mesh, bc);
        std::cout << "Problem size: " << mesh.nDof << " d.o.f.\n";

        //-----global mass matrix computation
        Eigen::SparseMatrix<Settings::precision> mG (mesh.nDof, mesh.nDof);
        mMat::computeGlobalM(mG, mesh);

        //-----Solving system of eqs. (full matrix method, no partitioning) 
        //     + Newton-Raphson + time integration (Newmark-beta method)

        ndofVectf force, dispOld, dispNew, velOld, velNew, accOld, accNew;
        
        Ic::setIC(dispOld, velOld, accOld, force, mG, mesh);
        bc.setBcInit(dispOld, force);
        dispNew = dispOld;
        velNew  = velOld;
        accNew  = accOld;

        Settings::precision time {Settings::t0 + Settings::dt};
        Settings::precision tSaveVTK {Settings::t0 + Settings::dtSaveVTK};
        int countVTK {0};
        Settings::precision normR {100.0};
        int newtonRaphsonCounter {0};

        clearPriorResults();

        constexpr const Settings::precision dt2 {Settings::dt*Settings::dt};
        constexpr const Settings::precision coeffM {1.0/(Settings::beta*dt2)};
        constexpr const Settings::precision coeffAcc {0.5-Settings::beta};

        ndofVectf forceInt = ndofVectf::Zero(mesh.nDof);
        ndofVectf residual;
        Eigen::SparseMatrix<Settings::precision> ktInt (mesh.nDof, mesh.nDof);
        FintKtInt::KtIntCache ktIntCache;
        Eigen::SparseMatrix<Settings::precision> ktEff (mesh.nDof, mesh.nDof);
        Settings::solverType solver;

        while (time < (Settings::tFinal + Settings::dt/2.0))
        {
            if (!Settings::importICFromStatic && (newtonRaphsonCounter == 0)){
                bc.setBcTimeDependent(force,time);
            }
            accNew = coeffM*(dispNew - dispOld - velOld*Settings::dt 
                - coeffAcc*dt2*accOld);
            velNew = velOld + (1.0-Settings::gamma)*Settings::dt*accOld
                + Settings::gamma*Settings::dt*accNew;

            FintKtInt::compute(mesh, dispNew, forceInt, ktInt, ktIntCache);

            residual.noalias() = mG*accNew;
            residual.noalias() += (forceInt - force);
            bc.embedBCintoResidual(residual);
            
            normR = residual.cwiseAbs().maxCoeff()/force.cwiseAbs().maxCoeff();
            
            if (normR < Settings::residualTargetTol){
                if (time > (tSaveVTK-Settings::dt/2.0)){
                    Export::exportResultsBinary(mesh, dispNew, countVTK, time);
                    tSaveVTK += Settings::dtSaveVTK;
                    ++countVTK;
                    std::cout << "Time instant " << time << " / " 
                        << Settings::tFinal << " s saved!\n";
                }
                std::swap(dispOld, dispNew);
                std::swap(velOld, velNew);
                std::swap(accOld, accNew);

                time += Settings::dt;
                newtonRaphsonCounter = 0;
                
                continue;
            }
            ktEff = coeffM * mG + ktInt;

            // as sparsity pattern of the ktEff matrix is fixed (assumption):
            if ((time < (Settings::t0 + 1.5*Settings::dt)) 
                && (newtonRaphsonCounter == 0))
            {
                solver.analyzePattern(ktEff);
                if(solver.info() != Eigen::Success) throw "Sparsity pattern analysis failed - linear solver.";
            }// evaluated only at the start of simulation-----------------------

            bc.embedBCintoKt(ktEff);

            solver.factorize(ktEff);
            if(solver.info()!=Eigen::Success) throw "Factorization failed - linear solver.";

            dispNew.noalias() -= solver.solve(residual);
            if (solver.info()!=Eigen::Success) throw "Linear solver failed.";

            if (++newtonRaphsonCounter > Settings::maxNewtonIters){
                throw "Maximum number of Newton-Raphson iterations reached!";
            }
        }

        std::cout << "\nResults saved in " << Settings::resultsDirName << '\n';
        std::cout << "To display time in paraview: filter -> python annotation"
            << " -> 'Time: %.3e' % (input.FieldData['TIME'][0])\n";
    }
    catch(const char* exceptionMsg){
        std::cerr << "Error: " << exceptionMsg << '\n';
    }
    catch(...){
        std::cerr << "Unrecognized error occured!\n";
    }

    return 0;
}