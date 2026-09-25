function Mel = elMassMat(coords,rho,t)

% computation of the element mass matrix using Gauss
% quadrature/integration.
% Inputs: - the nodal coordinates (x,y), 4x2 matrix
%         - rho, material density
%         - t, the thickness of the elements multiplying the integral

% full integration (which would give exact solution of the integral in the
% case in which the isoparamteric element coincides with the undistorted
% parent element -> Jacobian ratio = 1)

Mel = zeros(numel(coords)*[1,1]);
w = 1; % gaussian weights

for xi=[-1/sqrt(3), 1/sqrt(3)]
    for eta=[-1/sqrt(3), 1/sqrt(3)]
        Mel = Mel + w*w*(Nmat(xi,eta)' * rho * Nmat(xi,eta)).*detJ(coords,xi,eta);
    end
end

Mel = t*Mel;