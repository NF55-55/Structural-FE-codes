function Nmat = Nmatrix(xi,eta)

% Computation of the shape function matrix N
% Inputs: - xi, the value of the first isoparametric coordinate (Gauss/integr. point)
%         - eta, the value of the second isoparametric coordinate (Gauss/integr. point)

N1 = (1-xi)*(1-eta)/4;
N2 = (1+xi)*(1-eta)/4;
N3 = (1+xi)*(1+eta)/4;
N4 = (1-xi)*(1+eta)/4;

Nmat = [N1, 0, N2, 0, N3, 0, N4, 0;
        0, N1, 0, N2, 0, N3, 0, N4];

end

