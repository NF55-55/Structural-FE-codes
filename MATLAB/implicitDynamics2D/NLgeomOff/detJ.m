function jacobian = detJ(nodesCoord,xi,eta)

% computation of the determinant of the Jacobi's matrix, also called
% Jacobian.
% inputs: - nodes coordinates (x,y), 4x2 matrix
%         - xi, the value of the first isoparametric coordinate (Gauss/integr. point)
%         - eta, the value of the second isoparametric coordinate (Gauss/integr. point)

a = dxydxi(nodesCoord(:,1),eta);
b = dxydxi(nodesCoord(:,2),eta);
c = dxydeta(nodesCoord(:,1),xi);
d = dxydeta(nodesCoord(:,2),xi);

jacobian = (a*d - b*c);

end

