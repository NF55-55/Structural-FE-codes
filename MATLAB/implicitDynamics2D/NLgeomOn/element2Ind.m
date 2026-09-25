function indices = element2Ind(elNodes)

% Returns the indices of all dof of the element whose nodes are provided.
% The indices are returned as follows: [ind1x,ind1y,ind2x,ind2y,ind3x,...]
% following the order given by the connectivity table (CCW)
% Inputs:       - elNodes: IDs of element's nodes

indices = [2*elNodes(1)-1; 2*elNodes(1);
           2*elNodes(2)-1; 2*elNodes(2);
           2*elNodes(3)-1; 2*elNodes(3);
           2*elNodes(4)-1; 2*elNodes(4);];

end

