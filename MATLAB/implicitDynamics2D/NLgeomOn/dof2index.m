function nodesIndex = dof2index(nodeNumb,dir)

if strcmp(dir,'x')
    nodesIndex = 2*nodeNumb - 1;
elseif strcmp(dir,'y')
    nodesIndex = 2*nodeNumb;
else
    error("[error] Specified nodal direction is not 'x' or 'y'!")
end

end

