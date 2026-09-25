function Fint = FintCompute(S,nodes,elements,dk1,t)

% Computes the vector of internal nodal forces
% Inputs:   - S: material constitutive stiffness matrix
%           - nodes: initial nodal coordinates (the one at t0)
%           - elements: connectivity table (nodes ID per element)
%           - dk1: current (supposed) nodal displacement vector d_{k+1}
%           - t: elements thickness

Fint = zeros(2*size(nodes,1),1); % initializing the output vector
w = 1; % gaussian weights

% looping through every element of the mesh
for i = 1:size(elements,1)
    % current nodal disp. vector for element i
    dk1El = dk1(element2Ind(elements(i,:)));
    % original (at t0) nodal coordinate of element i (using total Lagrangian)
    nodesCoord0 = nodes(elements(i,:),:);
    % current nodal coordinates
    nodesCoord = nodes(elements(i,:),:)';
    nodesCoord = nodesCoord(:) + dk1El;
    nodesCoord = reshape(nodesCoord,2,[]);
    nodesCoord = nodesCoord';
    % looping through the gauss points of the i-th element
    for xi=[-1/sqrt(3), 1/sqrt(3)]
        for eta=[-1/sqrt(3), 1/sqrt(3)]
            invJ0 = invJac(nodesCoord0,xi,eta);
            % linear part of B_0
            B_L = [1,0,0,0; 0,0,0,1; 0,1,1,0]*invJ0*dN(xi,eta);
            % % E11 of non-linear part of B_0
            % B_NL1 = 0.5*dk1El'*dN(xi,eta)'*invJ0'*matE11*invJ0*dN(xi,eta);
            % % E22 of non-linear part of B_0 
            % B_NL2 = 0.5*dk1El'*dN(xi,eta)'*invJ0'*matE22*invJ0*dN(xi,eta);
            % % E12 of non-linear part of B_0
            % B_NL3 = 0.5*dk1El'*dN(xi,eta)'*invJ0'*matE12*invJ0*dN(xi,eta);
            % % Assemblying them
            % B_NL  = [B_NL1; B_NL2; B_NL3]; % 3x8
            % % Finding B_0
            % B_0 = B_L + B_NL;
            % Find stress vector P(u_{k+1})
            invJ  = invJac(nodesCoord,xi,eta);
            gradu = invJ*dN(xi,eta)*dk1El;
            E = [gradu(1) + 0.5*(gradu(1)^2 + gradu(3)^2);
                 gradu(4) + 0.5*(gradu(2)^2 + gradu(4)^2);
                 gradu(2) + gradu(3) + 0.5*(gradu(1)*gradu(2) + gradu(3)*gradu(4))];
            % Part of the gauss quadrature
            Fint(element2Ind(elements(i,:))) = Fint(element2Ind(elements(i,:))) + ...
                t*w* B_L' * S * E * detJ(nodesCoord,xi,eta);
        end
    end
end


end

