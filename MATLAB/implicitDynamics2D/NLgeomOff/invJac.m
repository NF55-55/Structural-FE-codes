function invJ = invJac(nodesCoord,xi,eta)

% computation of the inverse of the Jacobi's matrix for an element
% inputs: - nodes coordinates (x,y), 4x2 matrix
%         - xi, the value of the first isoparametric coordinate (Gauss/integr. point)
%         - eta, the value of the second isoparametric coordinate (Gauss/integr. point)

a = dxydxi(nodesCoord(:,1),eta);
b = dxydxi(nodesCoord(:,2),eta);
c = dxydeta(nodesCoord(:,1),xi);
d = dxydeta(nodesCoord(:,2),xi);

invJ = 1/(a*d-b*c)*[d,-b,0,0; -c,a,0,0; 0,0,d,-b; 0,0,-c,a];

end

