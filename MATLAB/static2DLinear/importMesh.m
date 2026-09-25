function [nodes,elements] = importMesh()

fprintf("\nNOTICE: For now, Abaqus and Ansys .inp files as well as SALOME .dat mesh\n"+...
    " files are supported, rename these so that they include the name of the program,\n"+...
    " written precisely as Abaqus, Ansys or SALOME (Example: 'modelAbaqus.inp').\n\n")

dirInp = ".\MESH\";
% check if the directory in which .inp should be stored exists
if ~exist(dirInp,"dir")
    mkdir(dirInp)
    error("[error] No '.\MESH' directory exist. Creating it... Now add the .inp file.")
end

inpFiles = dir(fullfile(dirInp,'*.inp'));
inpFiles2 = dir(fullfile(dirInp,'*.dat'));
inpFiles = struct2cell(inpFiles);
inpFiles2 = struct2cell(inpFiles2);
if isempty(inpFiles)
    if ~isempty(inpFiles2)
        inpFiles = {inpFiles2{1}};
    end
else
    if isempty(inpFiles2)
        inpFiles = {inpFiles{1}};
    else
        inpFiles = {inpFiles{1}, inpFiles2{1}};
    end
end
% check if inside the directory there are .MESH, if they are found see
% if it's an Abaqus or Ansys .inp file, based on the name given by the
% user, if both are found, Abaqus one is preferred.
if isempty(inpFiles)
    error("[error] No .MESH in the '.\MESH\' directory, please add the file.")
else
    nInpFiles = size(inpFiles,2);
    fileProgram = "";
    for i=1:nInpFiles
        inpFileName = inpFiles{1,i}; % first row is "name" category
        if contains(inpFileName,'Abaqus')
            fileProgram = 'Abaqus';
            fprintf("\nAbaqus file recognized in the './MESH/' directory!\n")
            break;
        elseif contains(inpFileName,'Ansys')
            fileProgram = 'Ansys';
            fprintf("\nAnsys file recognized in the './MESH/' directory!\n")
            break;
        elseif contains(inpFileName, 'SALOME')
            fileProgram = 'SALOME';
            fprintf("\nSALOME file recognized in the './MESH/' directory!\n")
        else
            if i == nInpFiles
                error("[error] The .inp/.dat file name doesn't contain neither" + ...
                    " 'Abaqus', 'Ansys' nor 'SALOME', rename the files if you " + ...
                    "forgot to do so!")
            end
        end
    end
end

% extract nodes coordinates and elements nodes from the .inp file based on
% the program (Ansys or Abaqus)
if strcmp(fileProgram,'Abaqus')
    [nodes, elements] = extractMeshFromAbaqusInp(dirInp,inpFileName);
    fprintf("\nAbaqus Mesh imported successfully!\n")
elseif strcmp(fileProgram,'Ansys')
    [nodes, elements] = extractMeshFromAnsysInp(dirInp,inpFileName);
    fprintf("\nAnsys Mesh imported successfully!\n")
else % SALOME is the only remaining alternative if code has reach this point
    [nodes, elements] = extractMeshFromSALOMEdat(dirInp,inpFileName);
    fprintf("\nSALOME Mesh imported successfully!\n")
end

end