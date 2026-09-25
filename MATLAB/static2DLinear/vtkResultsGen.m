function vtkResultsGen(nodes, elements, d, S)

% creating the .vtk file to visualize results (displacement and stress
% field) in Paraview

fprintf("\nWriting the results in .vtk file to visualize them in Paraview...\n")

numElem = size(elements, 1);
fid = fopen('.\VTK RESULTS\results.vtk', 'w');
fprintf(fid, ['# vtk DataFile Version 3.0\nFEA Unaveraged\nASCII\n' ...
    'DATASET UNSTRUCTURED_GRID\n\n']);

% 1. POINTS: Every element gets 4 unique points (Total = 4 * numElem), to
% be able to represent unaveraged stress field
fprintf(fid, 'POINTS %d float\n', numElem * 4);
reordElements = elements';
reordElements = reordElements(:);
coords = nodes(reordElements, :); % nNodesx2 matrix
fprintf(fid, '%f %f 0.0\n', coords'); 

% 2. CELLS
fprintf(fid, '\nCELLS %d %d\n', numElem, numElem * 5);
idEl = reshape(0:1:4*numElem-1,4,numElem);
fprintf(fid, '4 %d %d %d %d\n', idEl);

% 3. CELL_TYPES
fprintf(fid, '\nCELL_TYPES %d\n', numElem);
fprintf(fid, '%d\n', repmat(9, numElem, 1));

% 4. POINT_DATA (Displacements and Extrapolated Stresses)
fprintf(fid, '\nPOINT_DATA %d\n', numElem * 4);

% -- DISPLACEMENT --
fprintf(fid, 'VECTORS Displacement float\n');
dReshaped = reshape(d,2,[]);
dEl = dReshaped(:,reordElements);
fprintf(fid, '%f %f 0.0\n', dEl);

% -- UNAVERAGED STRESS TENSOR --
% fprintf(fid, '\nTENSORS Stress float\n');
% for e = 1:numElem
%     % Get 4 Gauss point stresses for this element
%     % Sigma is 3x4 (Sx, Sy, Sxy)' for all 4 nodes
%     dTemp = d(2*elements(e,:)+[-1;0]);
%     Sigma = elStresses(nodes(elements(e,:),:),S,reshape(dTemp,[],1),'all');
%     Sigma = extrapSigmaNodes(Sigma);
% 
%     for n = 1:4
%         sx = Sigma(1,n); sy = Sigma(2,n); sxy = Sigma(3,n);
%         % Write the 3x3 symmetric stress tensor as {sxx, sxy, sxz, syx, syy, syz, ...}
%         fprintf(fid, '%f %f 0.0 %f %f 0.0 0.0 0.0 0.0\n', sx, sxy, sxy, sy);
%     end
% end

% -- UNAVERAGED STRESS SCALAR COMPONENTS --

SigmaV = zeros(3,numElem*4); % rows are the three components S11,S22,S12

for e = 1:numElem
    dTemp = d(2*elements(e,:)+[-1;0]);
    Sigma = elStresses(nodes(elements(e,:),:),S,reshape(dTemp,[],1),'all'); % [3x4]
    Sigma = extrapSigmaNodes(Sigma);
    SigmaV(:,e*4-3:e*4) = Sigma;
end

fprintf(fid, '\nSCALARS S11 float\nLOOKUP_TABLE default\n');
fprintf(fid, '%f\n', SigmaV(1,:));

fprintf(fid, '\nSCALARS S22 float\nLOOKUP_TABLE default\n');
fprintf(fid, '%f\n', SigmaV(2,:));

fprintf(fid, '\nSCALARS S12 float\nLOOKUP_TABLE default\n');
fprintf(fid, '%f\n', SigmaV(3,:));

fprintf(fid, '\n');
fclose(fid);

fprintf("\nFile 'results.vtk' created in directory './VTK RESULTS/'.\n")

end

