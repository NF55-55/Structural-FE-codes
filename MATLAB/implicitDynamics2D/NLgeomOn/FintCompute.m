function Fint = FintCompute(S,nodes,elements,dk1,t)

% Computes the vector of internal nodal forces.
% Inputs:   - S: material constitutive stiffness matrix
%           - nodes: initial nodal coordinates (the one at t0)
%           - elements: connectivity table (nodes ID per element)
%           - dk1: current (supposed) nodal displacement vector d_{k+1}
%           - t: elements thickness

Fint = zeros(2*size(nodes,1),1); % Initializing Fint

w = 1; % gaussian weights

% looping through every element of the mesh
for i = 1:size(elements,1)
    % current nodal disp. vector for element i
    dk1El = dk1(element2Ind(elements(i,:)));
    % original (at t0) nodal coordinate of element i (using total Lagrangian)
    nodesCoord0 = nodes(elements(i,:),:)';
    % current nodal coordinates
    eulerCoord = nodesCoord0(:) + dk1El;
    % looping through the gauss points of the i-th element
    for xi=[-1/sqrt(3), 1/sqrt(3)]
        for eta=[-1/sqrt(3), 1/sqrt(3)]
            % deformation gradient tensor in material coordinates X_i
            [B0,F,detJac] = FandB0computeOptim(nodesCoord0',xi,eta,eulerCoord);
            % Green-Lagrange strain tensor
            E = 1/2 * (F' * F - eye(size(F)));
            % Green-Lagrange strain tensor (Voigt notation)
            E = [E(1,1); E(2,2); 2*E(1,2)];
            % Piola-Kirchhoff second stress tensor (Voigt notation)
            PK2 = S * E;
            % Piola-Kirchhoff second stress tensor
            PK2 = [PK2(1), PK2(3);
                   PK2(3), PK2(2)];
            % Nominal stress tensor
            P   = F * PK2;
            % Nominal stress in vectorial format: [Pxx; Pxy; Pyx; Pyy] see
            % notes (to be consistent with the derivative matrix)
            P = [P(1,1); P(1,2); P(2,1); P(2,2)];
            % Matrix 4 (nodes) x 2 (dir.)
            Fint(element2Ind(elements(i,:))) = Fint(element2Ind(elements(i,:))) + ...
                        (t * w*w * B0' * P * detJac);
        end
    end
end


end

