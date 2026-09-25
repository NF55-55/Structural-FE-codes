function ddeta = dxydeta(coord,xiVal)
% evaluated dx/deta or equivalently dy/deta
% inputs are: - x1,x2,x3,x4 or equivalently y1,y2,y3,y4 coordinates of the
%               four nodes of the element (CCW numbering of the nodes).
%             - the position xi in which the user wants to evaluate the
%               derivative, xi is the first isoparametric coordinate.
% For more information see the notes of the FEM derivation.

ddeta = 1/4*(-coord(1)-coord(2)+coord(3)+coord(4)+xiVal*(coord(1)-coord(2)+...
    coord(3)-coord(4)));

end

