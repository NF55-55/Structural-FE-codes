function dNmat = dN(xi,eta)

% Computation of the derivative of the shape function matrix w.r.t. the
% isoparametric coordinates xi and eta.
% Inputs: - xi, the value of the first isoparametric coordinate (Gauss/integr. point)
%         - eta, the value of the second isoparametric coordinate (Gauss/integr. point)

dNmat = (1/4)*[(-1+eta), 0, (1-eta), 0, (1+eta), 0, (-1-eta), 0;
               (-1+xi), 0, (-1-xi), 0, (1+xi), 0, (1-xi), 0;
               0, (-1+eta), 0, (1-eta), 0, (1+eta), 0,(-1-eta);
               0, (-1+xi), 0, (-1-xi), 0, (1+xi), 0, (1-xi)];

end

