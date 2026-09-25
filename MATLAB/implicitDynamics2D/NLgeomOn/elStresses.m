function SigmagaussPt = elStresses(coords0,S,eulerCoords)

% the output is the Cauchy stress sigma = {sigma_x; sigma_y; tau_xy}
% Inputs: - the nodal coordinates (x,y), 4x2 matrix
%         - S, the material constitutive stiffness matrix

xi = 1/sqrt(3)*[-1,1,1,-1];
eta = 1/sqrt(3)*[-1,-1,1,1];

SigmagaussPt = zeros(3,4);

for intPoint = 1:4
    % displacement gradient tensor in material coordinates X_i
    F = FcomputeOptim(coords0,xi(intPoint),eta(intPoint),eulerCoords);
    % Green-Lagrange strain tensor
    E = 1/2 * (F' * F - eye(size(F)));
    % Green-Lagrange strain tensor (Voigt notation)
    E = [E(1,1); E(2,2); 2*E(1,2)];
    % Piola-Kirchhoff second stress tensor (Voigt notation)
    %PK2gaussPt(:,intPoint) = S * E;
    PK2gaussPt = S * E;
    PK2gaussPt = [PK2gaussPt(1), PK2gaussPt(3);
                  PK2gaussPt(3), PK2gaussPt(2)];
    sigmaTemp = (1/det(F)) * F * PK2gaussPt * F';
    SigmagaussPt(:,intPoint) = [sigmaTemp(1,1); sigmaTemp(2,2); sigmaTemp(1,2)];
end

end

