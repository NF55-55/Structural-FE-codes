function ddxi = dxydxi(coord,etaVal)
% evaluated dx/dxi or equivalently dy/dxi
% inputs are: - x1,x2,x3,x4 or equivalently y1,y2,y3,y4 coordinates of the
%               four nodes of the element (CCW numbering of the nodes).
%             - the position eta in which the user wants to evaluate the
%               derivative, eta is the second isoparametric coordinate.
% For more information see the notes of the FEM derivation.

ddxi = 1/4*(-coord(1)+coord(2)+coord(3)-coord(4)+etaVal*(coord(1)-coord(2)+...
    coord(3)-coord(4)));

end

