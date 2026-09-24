#include "export.h"

namespace Export{
    void exportResultsASCII(const Mesh& mesh, const DispEigen& disp)
    {
        std::cout << "Writing ./'VTK RESULTS'/results.vtk to visualize the "
            << "output fields (UNAVERAGED) in Paraview...";
        std::ofstream resultsFile {"./VTK RESULTS/results.vtk"};
        if (!resultsFile){
            throw "Unable to open ./'VTK RESULTS'/results.vtk - export.cpp";
        }
        std::ios_base::sync_with_stdio(false); // de-sync with C streams

        resultsFile << "# vtk DataFile Version 3.0\nFEA Unaveraged\nASCII\n"
            << "DATASET UNSTRUCTURED_GRID\n\n";
        resultsFile << "POINTS " << mesh.nElements*4 << " float\n";
        Eigen::IOFormat coordsFormat;
        if (Settings::fullPrecisionASCIIWriting){
            coordsFormat = {Eigen::FullPrecision, Eigen::DontAlignCols, " "
                , "\n", "", " 0.0"};
        } else{
            coordsFormat = {Eigen::StreamPrecision, Eigen::DontAlignCols, " "
                , "\n", "", " 0.0"};
        }
        Eigen::Vector<int,Eigen::Dynamic> flattenedNodesIDEl = mesh.elements
            .reshaped<Eigen::RowMajor>();
        resultsFile << (mesh.nodes(flattenedNodesIDEl,Eigen::placeholders::all))
            .format(coordsFormat) << '\n';

        resultsFile << "\nCELLS " << mesh.nElements << ' ' << mesh.nElements*5 
            << "\n";
        Eigen::IOFormat cellsNodesFormat {Eigen::StreamPrecision
            ,Eigen::DontAlignCols, " ", "\n", "4 "};
        Eigen::Vector<int,Eigen::Dynamic> tempElIndx = Eigen::Vector<int
            ,Eigen::Dynamic>::LinSpaced(mesh.nElements*Settings::nNodesEl, 0
            , mesh.nElements*Settings::nNodesEl-1);
        Eigen::Matrix<int,Eigen::Dynamic,Settings::nNodesEl> reshapedTempElIndx 
            = tempElIndx.reshaped<Eigen::RowMajor>(Eigen::AutoSize,Settings::
            nNodesEl);
        resultsFile << reshapedTempElIndx.format(cellsNodesFormat) << '\n';

        resultsFile << "\nCELL_TYPES " << mesh.nElements << '\n';
        Eigen::IOFormat cellTypeFormat {Eigen::StreamPrecision
            , Eigen::DontAlignCols, "", "\n"};
        Eigen::Vector<int,Eigen::Dynamic> tempCellType = Eigen::Vector<int
            ,Eigen::Dynamic>::Constant(mesh.nElements, 9);
        resultsFile << tempCellType.format(cellTypeFormat) << '\n';

        resultsFile << "\nPOINT_DATA " << mesh.nElements*4 << '\n';
        resultsFile << "VECTORS Displacement float\n";
        Eigen::IOFormat dispFormat;
        if (Settings::fullPrecisionASCIIWriting){
            dispFormat = {Eigen::FullPrecision, Eigen::DontAlignCols, " "
                , "\n", "", " 0.0"};
        } else{
            dispFormat = {Eigen::StreamPrecision, Eigen::DontAlignCols, " ", "\n"
                , "", " 0.0"};
        }
        Eigen::Matrix<Settings::precision,Eigen::Dynamic,2> reshapedDisp 
            = disp.reshaped<Eigen::RowMajor>(Eigen::AutoSize,2);
        resultsFile << reshapedDisp(flattenedNodesIDEl,Eigen::placeholders::all)
            .format(dispFormat) << '\n';

        SigmaMat sigmaM (4*mesh.nElements, stressComp);
        computeStressField(mesh, reshapedDisp, sigmaM);
        Eigen::IOFormat stressFormat;
        if (Settings::fullPrecisionASCIIWriting){
            stressFormat = {Eigen::FullPrecision, Eigen::DontAlignCols, ""
                , "\n", "", ""};
        } else{
            stressFormat = {Eigen::StreamPrecision, Eigen::DontAlignCols, ""
                , "\n", "", ""};
        }
        resultsFile << "\nSCALARS S11 float\nLOOKUP_TABLE default\n";
        resultsFile << sigmaM(Eigen::placeholders::all,0).format(stressFormat)
            << '\n';
        resultsFile << "\nSCALARS S22 float\nLOOKUP_TABLE default\n";
        resultsFile << sigmaM(Eigen::placeholders::all,1).format(stressFormat)
            << '\n';
        resultsFile << "\nSCALARS S12 float\nLOOKUP_TABLE default\n";
        resultsFile << sigmaM(Eigen::placeholders::all,2).format(stressFormat)
            << '\n';

        resultsFile.close();
        std::cout << "Done!\n\n";
    }

    void exportResultsBinary(const Mesh& mesh, const DispEigen& disp){
        std::cout << "Writing ./'VTK RESULTS'/results.bin.vtk to visualize the "
            << "output fields (UNAVERAGED) in Paraview...\n";
        std::ofstream file("./VTK RESULTS/results.bin.vtk", std::ios::binary);
        if (!file) throw std::runtime_error("Unable to open VTK output file.");

        // Increasing file buffer to 1MB to reduce OS writing calls
        std::vector<char> buf(1024 * 1024);
        file.rdbuf()->pubsetbuf(buf.data(), buf.size());

        file << "# vtk DataFile Version 3.0\nFEA Unaveraged\nBINARY\nDATASET UNSTRUCTURED_GRID\n\n";

        const int nPoints = mesh.nElements * Settings::nNodesEl;

        // string used to declare the floating point precision
        std::string_view floatingTypeStr = (isDouble(Settings::sqrt3) 
            ? "double" : "float");

        file << "POINTS " << nPoints << ' ' << floatingTypeStr << '\n';
        std::vector<Settings::precision> pointData(nPoints * 3);
        
        int idx = 0;
        for (int e=0; e<mesh.nElements; ++e) {
            for (int n=0; n<static_cast<int>(Settings::nNodesEl); ++n) {
                int nodei = mesh.elements(e, n);
                pointData[idx++] = swapEndian(mesh.nodes(nodei, 0));
                pointData[idx++] = swapEndian(mesh.nodes(nodei, 1));
                pointData[idx++] = swapEndian(0.0);
            }
        }
        file.write(reinterpret_cast<const char*>(pointData.data())
            , pointData.size() * sizeof(Settings::precision));
        // In the creation of binary vtk files, data is passed from memory as a 
        // contiguous data container (efficient), as the total size and the 
        // size of each set is specified (or assumed as in this case, 3D space
        // coordiantes = size 3).
        file << "\n";

        int CellEntries {mesh.nElements * static_cast<int>(Settings::nNodesEl 
            + 1uz)};
        file << "CELLS " << mesh.nElements << " " << CellEntries << "\n";
        std::vector<int32_t> cellData(CellEntries);
        idx = 0;
        int pointCounter = 0;
        for (int e=0; e<mesh.nElements; ++e) {
            cellData[idx++] = swapEndian<int32_t>(Settings::nNodesEl);
            for (int n=0; n<static_cast<int>(Settings::nNodesEl); ++n) {
                cellData[idx++] = swapEndian<int32_t>(pointCounter++);
            }
        }
        file.write(reinterpret_cast<const char*>(cellData.data())
            , cellData.size() * sizeof(int32_t));
        file << "\n";

        file << "CELL_TYPES " << mesh.nElements << "\n";
        // VTK 9 is linear (4-node) quad element
        std::vector<int32_t> cellTypes(mesh.nElements, swapEndian<int32_t>(9));
        file.write(reinterpret_cast<const char*>(cellTypes.data())
            , cellTypes.size() * sizeof(int32_t));
        file << "\n";

        file << "POINT_DATA " << nPoints << "\n";
        file << "VECTORS Displacement " << floatingTypeStr << '\n';
        std::vector<Settings::precision> dispData(nPoints * 3);
        idx = 0;
        Eigen::Matrix<Settings::precision, Eigen::Dynamic, Settings::dofNode> 
            reshapedDisp = disp.reshaped<Eigen::RowMajor>(Eigen::AutoSize
            , Settings::dofNode);
        
        for (int e=0; e<mesh.nElements; ++e) {
            for (int n=0; n<static_cast<int>(Settings::nNodesEl); ++n) {
                int nodeID = mesh.elements(e, n);
                dispData[idx++] = swapEndian(reshapedDisp(nodeID, 0));
                dispData[idx++] = swapEndian(reshapedDisp(nodeID, 1));
                dispData[idx++] = swapEndian(0.0f);
            }
        }
        file.write(reinterpret_cast<const char*>(dispData.data())
            , dispData.size() * sizeof(Settings::precision));
        file << "\n";

        SigmaMat sigmaM(4 * mesh.nElements, stressComp);
        computeStressField(mesh, reshapedDisp, sigmaM);
        // lambda function to write stress components
        auto writeStressScalar = [&](const char* name, int col){
            file << "SCALARS " << name << ' ' << floatingTypeStr << '\n' 
                << "LOOKUP_TABLE default\n";
            std::vector<Settings::precision> stressData(nPoints);
            for (int i=0; i<nPoints; ++i){
                stressData[i] = swapEndian(sigmaM(i, col));
            }
            file.write(reinterpret_cast<const char*>(stressData.data())
                , stressData.size() * sizeof(Settings::precision));
            file << "\n";
        };
        writeStressScalar("S11", 0);
        writeStressScalar("S22", 1);
        writeStressScalar("S12", 2);
    }

    void computeStressField(const Mesh& mesh, DispMatEigen& dispMat
        , SigmaMat& sigmaM)
    {
        Settings::precision dxdxi   {};
        Settings::precision dydxi   {};
        Settings::precision dxdeta  {};
        Settings::precision dydeta  {};
        Settings::precision detJ    {};
        Settings::precision invDetJ {};
        Settings::precision xi {};
        Settings::precision eta {};
        Settings::precision r {};
        Settings::precision s {};
        // extrapolation basis functions vector
        Eigen::Vector<Settings::precision,Settings::nNodesEl> Nextr;
        NodesCoordEigen coords;
        dNEigen dN;
        InvJEigen invJ;
        BEigen bEl;
        Eigen::Vector<Settings::precision,Settings::nDofEl> dispEl;
        Eigen::Matrix<Settings::precision,Settings::nNodesEl,stressComp> 
            sigmaTemp;
        for (int e{0}; e<mesh.nElements; ++e){
            dispEl.noalias() = dispMat(mesh.elements(e
                ,Eigen::placeholders::all),Eigen::placeholders::all)
                .reshaped<Eigen::RowMajor>();
            coords = mesh.nodes(mesh.elements.row(e),Eigen::placeholders::all);
            for (int intPtCCW{0}; intPtCCW<Settings::nNodesEl; ++intPtCCW){
                xi = xiV[intPtCCW];
                eta = etaV[intPtCCW];
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
                sigmaTemp(intPtCCW,Eigen::placeholders::all).noalias() 
                    = Material::C * bEl * dispEl;
            }
            // Extrapolating stresses at the nodes starting from the one found 
            // for the integration points
            for (int nodei{0}; nodei<Settings::nNodesEl; ++nodei){
                r = rV[nodei];
                s = sV[nodei];
                Nextr(0) = (1-r)*(1-s)/4.0;
                Nextr(1) = (1+r)*(1-s)/4.0;
                Nextr(2) = (1+r)*(1+s)/4.0;
                Nextr(3) = (1-r)*(1+s)/4.0;
                sigmaM(e*Settings::nNodesEl+nodei,Eigen::placeholders::all)
                    .noalias() = (Nextr.asDiagonal() * sigmaTemp).colwise()
                    .sum();
            }
        }
    }
}