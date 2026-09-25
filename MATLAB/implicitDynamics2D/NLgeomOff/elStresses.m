function sigma = elStresses(coords,S,dEl,typeSigma)

% the output is sigma = {sigma_x; sigma_y; tau_xy}
% Inputs: - the nodal coordinates (x,y), 4x2 matrix
%         - S, the material constitutive stiffness matrix

sigma = zeros(3,4); % 4 given the number of Gauss points in this linear el.
xi = 1/sqrt(3)*[-1,1,1,-1];
eta = 1/sqrt(3)*[-1,-1,1,1];

for intPoint = 1:4
    Bel = [1,0,0,0; 0,0,0,1; 0,1,1,0] * invJac(coords,xi(intPoint),eta(intPoint)) * ...
        dN(xi(intPoint),eta(intPoint));
    sigma(:,intPoint) = S * Bel * dEl;
end

if strcmp(typeSigma,'x')
    sigma = sigma(1,:);
elseif strcmp(typeSigma,'y')
    sigma = sigma(2,:);
elseif strcmp(typeSigma,'xy')
    sigma = sigma(3,:);
elseif strcmp(typeSigma,'all')
    % just leave sigma full as it is
else
    error("[error] The type of stress requested is invalid. Valid ones: 'x','y','xy','all'")
end

end

