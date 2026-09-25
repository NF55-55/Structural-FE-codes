function [B0,F,detJac] = FandB0computeOptim(coord,xi,eta,ec)

% returning both B0, det(J) and F as to compute F we also need B0
% improvement of B0computeOptim and FcomputeOptim
% Inputs:   - coords, the Lagrangian (in undeformed state) coordinates 
%             of the element's nodes (4x2)
%           - xi + eta, the parent coordinates of isoparametric formulation
%           - ec, the Eulerian (current deformed state) coordinates
%             of the element's nodes (x = X + u_x) (8x1) 

a = 1/4*(-coord(1,1)+coord(2,1)+coord(3,1)-coord(4,1)+eta*(coord(1,1)-coord(2,1)+...
    coord(3,1)-coord(4,1)));
b = 1/4*(-coord(1,2)+coord(2,2)+coord(3,2)-coord(4,2)+eta*(coord(1,2)-coord(2,2)+...
    coord(3,2)-coord(4,2)));
c = 1/4*(-coord(1,1)-coord(2,1)+coord(3,1)+coord(4,1)+xi*(coord(1,1)-coord(2,1)+...
    coord(3,1)-coord(4,1)));
d = 1/4*(-coord(1,2)-coord(2,2)+coord(3,2)+coord(4,2)+xi*(coord(1,2)-coord(2,2)+...
    coord(3,2)-coord(4,2)));

detJac = (a*d-b*c);

invJ = 1/detJac*[d,-b,0,0; -c,a,0,0; 0,0,d,-b; 0,0,-c,a];
dN   = (1/4)*[(-1+eta), 0, (1-eta), 0, (1+eta), 0, (-1-eta), 0;
              (-1+xi), 0, (-1-xi), 0, (1+xi), 0, (1-xi), 0;
              0, (-1+eta), 0, (1-eta), 0, (1+eta), 0, (-1-eta);
              0, (-1+xi), 0, (-1-xi), 0, (1+xi), 0, (1-xi)];

B0 = invJ * dN;

F = B0 * ec; F = [F(1),F(2); F(3),F(4)];

end

