function convBCandF(datF,nameBCorF,fullPath,targetPath)

if sum(strcmp(datF,nameBCorF+".dat"))>=1
    file = dir(fullPath+nameBCorF+".dat");
    if file.bytes ~= 0
        lenLines = readmatrix(fullPath+nameBCorF+".dat",'Delimiter', ' ', ...
            'ConsecutiveDelimitersRule', 'join', 'LeadingDelimitersRule', 'ignore', ...
            'TrailingDelimitersRule', 'ignore','Range',[1,1,1,1],'FileType','text') + 1;
        nodeId = readmatrix(fullPath+nameBCorF+".dat",'Delimiter', ' ', ...
            'ConsecutiveDelimitersRule', 'join', 'LeadingDelimitersRule', 'ignore', ...
            'TrailingDelimitersRule', 'ignore','Range',[2,1,lenLines,1], 'FileType','text');
        fid = fopen(targetPath+nameBCorF+".txt","w");
        if lenLines == 2
            fprintf(fid,'%d',nodeId(1));
        else
            fprintf(fid,'%d, ',nodeId(1:end-1));
            fprintf(fid,'%d',nodeId(end));
        end
        fclose(fid);
    else
        fid = fopen(targetPath+nameBCorF+".txt","w");
        fprintf(fid,'');
        fclose(fid);
    end
else
    strErr = "[error] No file named"+nameBCorF+", check names.";
    error(strErr)
end

end

