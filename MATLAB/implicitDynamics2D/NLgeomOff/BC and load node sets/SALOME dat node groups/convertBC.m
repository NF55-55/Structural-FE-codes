targetPath = "..\";
fullPath = ".\";
datFiles = dir(fullfile(fullPath,"*.dat"));

if size(datFiles,1)~=4
    error("[error] All BC and load node group files should be stored in './BC and " + ...
        "load node set/SALOME dat node groups/[BCx BCy Fx Fy].dat', even if empty.")
end
datFiles = struct2cell(datFiles);
datF = strings(4,1);
for i = 1:4
    datF(i) = datFiles{1,i};
end

% converting .dat to .txt with format used in simulation
fprintf("\nMake sure that the .dat files with no data have no white spaces!\n")
nameBCorF = ["BCx","BCy","Fx","Fy"];
for i = 1:4
    convBCandF(datF,nameBCorF(i),fullPath,targetPath)
end