function [nodes, elements] = extractMeshFromSALOMEdat(dirInp,inpFileName)

lines = char(readlines(dirInp+inpFileName));
% reading the number of nodes and elements reported in first line
nNodesEl = readmatrix(dirInp+inpFileName,'Delimiter',' ','Range', [1,1,1,2]); 
bottomFile = sum(nNodesEl) + 1; % adding 1 row as the first one is for mesh info
headers = find(sum(lines(:,1:2)==['1' ' '],2)==2);
nodeHeaderLine = headers(1);
nodeBottomLine = headers(2)-1;
clear lines

% extracting nodal coordinates
nodes = readmatrix(dirInp+inpFileName,'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join', 'LeadingDelimitersRule', 'ignore', ...
    'TrailingDelimitersRule', 'ignore','Range',[nodeHeaderLine,2,nodeBottomLine,3], ...
    'FileType','text');

% extracting elements' nodes number
elementsType = readmatrix(dirInp+inpFileName, 'Delimiter', ' ', 'Range', ...
    [nodeBottomLine+1,2,bottomFile,2]);
% 204 is the ID for quad linear element
elemHeaderLine = find(elementsType==204,1,"first") + nodeBottomLine;
elemBottomLine = find(elementsType==204,1,"last") + nodeBottomLine;
elements = readmatrix(dirInp+inpFileName,'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join', 'LeadingDelimitersRule', 'ignore', ...
    'TrailingDelimitersRule', 'ignore','Range',[elemHeaderLine,3,elemBottomLine,6], ...
    'FileType','text');

% Finally the group of nodes saved in '.\BC and load node sets\SALOME dat node groups\'
% are converted into the actualy BC and load nodes files used in
% simulation:
run(".\BC and load node sets\SALOME dat node groups\convertBC.m")

end

