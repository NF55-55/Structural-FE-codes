function [nodes, elements] = extractMeshFromAnsysInp(dirInp,inpFileName)

lines = readlines(dirInp+inpFileName);
nodeHeaderLine = find(contains(lines,'*********** Nodes'),1);
nodeBottomLine = find(contains(lines,'*********** Elements'),1);
elementHeaderLine = nodeBottomLine;
elementBottomLine = find(contains(lines,' !  done creating elements'),1);
clear lines

% extracting nodal coordinates
nodes = readmatrix(dirInp+inpFileName,'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join', 'LeadingDelimitersRule', 'ignore', ...
    'TrailingDelimitersRule', 'ignore','Range',[nodeHeaderLine+3,2,nodeBottomLine-3,3], ...
    'FileType','text');

% extracting elements' nodes number
elements = readmatrix(dirInp+inpFileName,'Delimiter', ' ', ...
    'ConsecutiveDelimitersRule', 'join', 'LeadingDelimitersRule', 'ignore', ...
    'TrailingDelimitersRule', 'ignore','Range', ...
    [elementHeaderLine+6,2,elementBottomLine-4,5], 'FileType','text');

end

