function sigma = extrapSigmaNodes(sigma)

% extrapolating sigma at the nodes starting from its values at the Gauss
% points

% CCW 
r=sqrt(3)*[-1,1,1,-1];
s =sqrt(3)*[-1,-1,1,1];
sTemp = sigma;

for i = 1:4
    N1 = (1-r(i))*(1-s(i))/4;
    N2 = (1+r(i))*(1-s(i))/4;
    N3 = (1+r(i))*(1+s(i))/4;
    N4 = (1-r(i))*(1+s(i))/4;
    sigma(:,i) = N1*sTemp(:,1) + N2*sTemp(:,2) + N3*sTemp(:,3) + N4*sTemp(:,4);
end

end

