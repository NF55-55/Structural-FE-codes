function invJ = invJac(coord,xi,eta)

% computation of the inverse of the Jacobi's matrix for an element
% inputs: - nodes coordinates (x,y), 4x2 matrix
%         - xi, the value of the first isoparametric coordinate (Gauss/integr. point)
%         - eta, the value of the second isoparametric coordinate (Gauss/integr. point)

a = 1/4*(-coord(1,1)+coord(2,1)+coord(3,1)-coord(4,1)+eta*(coord(1,1)-coord(2,1)+...
    coord(3,1)-coord(4,1)));
b = 1/4*(-coord(1,2)+coord(2,2)+coord(3,2)-coord(4,2)+eta*(coord(1,2)-coord(2,2)+...
    coord(3,2)-coord(4,2)));
c = 1/4*(-coord(1,1)-coord(2,1)+coord(3,1)+coord(4,1)+xi*(coord(1,1)-coord(2,1)+...
    coord(3,1)-coord(4,1)));
d = 1/4*(-coord(1,2)-coord(2,2)+coord(3,2)+coord(4,2)+xi*(coord(1,2)-coord(2,2)+...
    coord(3,2)-coord(4,2)));

invJ = 1/(a*d-b*c)*[d,-b,0,0; -c,a,0,0; 0,0,d,-b; 0,0,-c,a];

end

