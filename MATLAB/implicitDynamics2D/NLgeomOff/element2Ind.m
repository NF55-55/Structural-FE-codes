function indices = element2Ind(elNodes)

% Returns the indices of all dof of the element whose nodes are provided.
% The indices are returned as follows: [ind1x,ind1y,ind2x,ind2y,ind3x,...]
% following the order given by the connectivity table (CCW)
% Inputs:       - elNodes: IDs of element's nodes

indices = [dof2index(elNodes,'x'); dof2index(elNodes,'y')];
indices = indices(:);

end

