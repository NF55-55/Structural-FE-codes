function [nodes, elements] = extractMeshFromAbaqusInp(dirInp,inpFileName)

lines = readlines(dirInp+inpFileName);
nodeHeaderLine = find(contains(lines,'*Node'),1);
elementHeaderLine = find(contains(lines,'*Element'),1);

% extracting nodal coordinates
strRange = num2str(int64(nodeHeaderLine+1))+":"+num2str(int64(elementHeaderLine-1));
nodesBlock = readmatrix(dirInp+inpFileName,'Delimiter',',','Range',strRange, ...
    'FileType','text');
nodes = nodesBlock(:,2:end);

% k = 1;
% for i=nodeHeaderLine+1:elementHeaderLine-1
%     tempStr = char(lines(i));
%     commaIndex = strfind(tempStr,',');
%     nodes(k,:) = [str2double(tempStr(commaIndex(1)+1:commaIndex(2)-1)), ...
%         str2double(tempStr(commaIndex(2)+1:end))];
%     k = k + 1;
% end

% extracting elements' nodes number
elementBottomLine = find(contains(lines(elementHeaderLine+1:end),'*'),1) + ...
     elementHeaderLine;

strRange = num2str(int64(elementHeaderLine+1))+":"+num2str(int64(elementBottomLine-1));
elemBlock = readmatrix(dirInp+inpFileName,'Delimiter',',','Range',strRange, ...
    'FileType','text');
elements = elemBlock(:,2:end);

% k = 1;
% for i=elementHeaderLine+1:elementBottomLine-1
%     tempStr = lines(i);
%     strVals = extract(tempStr, digitsPattern);
%     elements(k,:) = str2double(strVals(2:end));
%     k = k + 1;
% end

end

