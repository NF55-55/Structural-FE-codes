function [F0, KtInt] = KtIntNum(S,nodes,elements,dk1,t)

F0 = FintCompute(S,nodes,elements,dk1,t); % reference value for this time-step and 
% increment
delta = sqrt(eps); % variation in the nodal displacement for numerical computation of
% numerical internal tangent stiffness matrix initialization
KtInt = zeros(length(dk1),length(dk1));
for i = 1:length(nodes)*2
    incr = zeros(size(F0));
    incr(i) = delta;
    newdk1 = dk1 + incr;
    F1 = FintCompute(S,nodes,elements,newdk1,t);
    dFintdui = (F1-F0)/delta;
    KtInt(:,i) = dFintdui;
end

end

