function Kel = elStiffMat(coords,S,t)

% computation of the element stiffness matrix using Gauss
% quadrature/integration.
% Inputs: - the nodal coordinates (x,y), 4x2 matrix
%         - S, the material constitutive stiffness matrix
%         - t, the thickness of the elements multiplying the integral

% full integration (which would give exact solution of the integral in the
% case in which the isoparamteric element coincides with the undistorted
% parent element -> Jacobian ratio = 1)

Kel = zeros(numel(coords)*[1,1]);
w = 1; % gaussian weights

for xi=[-1/sqrt(3), 1/sqrt(3)]
    for eta=[-1/sqrt(3), 1/sqrt(3)]
        % element's strain-displacement matrix
        Bel = [1,0,0,0; 0,0,0,1; 0,1,1,0] * invJac(coords,xi,eta) * dN(xi,eta);
        Kel = Kel + t*w*w*(Bel' * S * Bel).*detJ(coords,xi,eta);
    end
end

end