#include <iostream>
#include <fstream>
#include <algorithm>
#include "Eigen/Eigen"
#include "Eigen/PardisoSupport"
#include "Mesh.h"
#include "import.h"
#include "settings.h"
#include "TimerHighRes.h"
#include "material.h"
#include "Bc.h"
#include "rearrangeDof.h"
#include "kmMatrices.h"
#include "export.h"
#include "Ic.h"
#include "clearResults.h"

int main(){
    //Timer timer {};
    //timer.start();
    Mesh mesh {};
    Bc bc {}; // boundary conditions class
    using ndofVectf = Eigen::Vector<Settings::precision, Eigen::Dynamic>;
    using ndofVecti = Eigen::Vector<int, Eigen::Dynamic>;
    using SizT = std::size_t;
    ndofVectf dispOld;
    ndofVectf force;
    
    // To handle exceptions in any function/method
    try{
        //-----Importing mesh and boundary conditions
        importMesh(mesh, bc);
        std::cout << "Problem size: " << mesh.nDof << " d.o.f.\n";

        //-----(global mass matrix + global stiffness matrix) computation
        Eigen::SparseMatrix<Settings::precision> mG (mesh.nDof, mesh.nDof);
        Eigen::SparseMatrix<Settings::precision> kG (mesh.nDof, mesh.nDof);
        kmMat::computekGmG(mesh, mG, kG);

        //-----Solving system of eqs. (partitioning method) 
        //     & time integration (Newmark-beta method)

        dispOld = ndofVectf::Constant(mesh.nDof,Settings::NaN);
        force = ndofVectf::Zero(mesh.nDof);

        bc.setBcInit(dispOld, force);

        // Re-arranging vectors (partitioning method) + saving mapping vector
        ndofVecti mapDof;
        mapDof = ndofVecti::LinSpaced(mesh.nDof, 0, mesh.nDof-1);
        rearrange(dispOld, force, mapDof);
        Eigen::PermutationMatrix<Eigen::Dynamic, Eigen::Dynamic, int> 
            mapDofMat {mapDof};
        kG = kG.twistedBy(mapDofMat.transpose()); // Re-arranging
        mG = mG.twistedBy(mapDofMat.transpose()); // Re-arranging

        int lastNanIndx {mesh.nDof-1};
        for (; lastNanIndx>0; --lastNanIndx){
            if (std::isnan(dispOld(lastNanIndx))) break;
        }
        int nUnknownDisp {lastNanIndx+1};
        int nKnownDisp {mesh.nDof - nUnknownDisp};

        ndofVectf dispNew, velOld, velNew, accOld, accNew, dispReord;
        
        Ic::setIC(dispOld, velOld, accOld, force, mG, kG, nUnknownDisp, mapDof);
        dispNew = dispOld;
        velNew  = velOld;
        accNew  = accOld;

        Settings::precision time {Settings::t0 + Settings::dt};
        Settings::precision tSaveVTK {Settings::t0 + Settings::dtSaveVTK};
        int countVTK {0};

        clearPriorResults();

        constexpr const Settings::precision dt2 {Settings::dt*Settings::dt};
        constexpr const Settings::precision coeffM {1.0/(Settings::beta*dt2)};
        constexpr const Settings::precision coeffAcc {0.5-Settings::beta};
        
        Eigen::SparseMatrix<Settings::precision> effStiffMat 
            (nUnknownDisp,nUnknownDisp);
        effStiffMat = coeffM*mG.block(0,0,nUnknownDisp,nUnknownDisp) 
            + kG.block(0,0,nUnknownDisp,nUnknownDisp);
        Settings::solverType solver;
        solver.compute(effStiffMat);
        if(solver.info()!=Eigen::Success) {
            throw "Linear solver failed decomposition.";
        }

        Eigen::Vector<Settings::precision,Eigen::Dynamic> rhsEquation;

        while (time < (Settings::tFinal + Settings::dt/2.0)){
            if (!Settings::importICFromStatic){
                bc.setBcTimeDependent(force,mapDof,time,nUnknownDisp);
            }

            rhsEquation.noalias() = force.head(nUnknownDisp);
            rhsEquation.noalias() += coeffM*mG.topLeftCorner(nUnknownDisp
                ,nUnknownDisp)*(dispOld.head(nUnknownDisp) 
                + Settings::dt*velOld.head(nUnknownDisp) 
                + coeffAcc*dt2*accOld.head(nUnknownDisp));
            rhsEquation.noalias() += coeffM*mG.topRightCorner(nUnknownDisp
                ,nKnownDisp)*(-dispNew.tail(nKnownDisp) 
                + dispOld.tail(nKnownDisp) 
                + Settings::dt*velOld.tail(nKnownDisp) 
                + coeffAcc*dt2*accOld.tail(nKnownDisp));
            rhsEquation.noalias() -= kG.topRightCorner(nUnknownDisp,nKnownDisp)
                *dispNew.tail(nKnownDisp);
            
            dispNew.head(nUnknownDisp) = solver.solve(rhsEquation);
            if(solver.info()!=Eigen::Success){
                throw "Linear solver failed solving(1).";
            }
            
            accNew.head(nUnknownDisp).noalias() = coeffM
                *(dispNew.head(nUnknownDisp) - dispOld.head(nUnknownDisp) 
                - velOld.head(nUnknownDisp)*Settings::dt - coeffAcc
                *dt2*accOld.head(nUnknownDisp));

            velNew.head(nUnknownDisp).noalias() = velOld.head(nUnknownDisp)
                + (1.0-Settings::gamma)*Settings::dt*accOld.head(nUnknownDisp)
                + Settings::gamma*Settings::dt*accNew.head(nUnknownDisp);

            if (time > (tSaveVTK-Settings::dt/2.0)){
                dispReord = mapDofMat * dispNew;
                Export::exportResultsBinary(mesh, dispReord, countVTK, time);
                tSaveVTK += Settings::dtSaveVTK;
                ++countVTK;
                std::cout << "Time instant " << time << " / " 
                    << Settings::tFinal << " s saved!\n";
            }

            std::swap(dispOld, dispNew);
            std::swap(velOld, velNew);
            std::swap(accOld, accNew);

            time += Settings::dt;
        }

        std::cout << "\nResults saved in " << Settings::resultsDirName << '\n';
        std::cout << "To display time in paraview: filter -> python annotation"
            << " -> 'Time: %.3e' % (input.FieldData['TIME'][0])\n";
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
