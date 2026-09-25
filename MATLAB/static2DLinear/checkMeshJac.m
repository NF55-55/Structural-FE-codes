function checkMeshJac(nodes,elements,plotBool)

nodesCoord = zeros(size(elements,2),size(nodes,2));
nodalJac = zeros(size(elements,1),size(elements,2)); % jacobian at the nodes initialization

for i=1:size(elements,1) % for each element doing the following
    nodesCoord = nodes(elements(i,:),:);
    k = 1;
    for xi=[-1,1]
        for eta=[-1,1]
            % for the following nomenclature, explanation is in the notes
            a = dxydxi(nodesCoord(:,1),eta);
            b = dxydxi(nodesCoord(:,2),eta);
            c = dxydeta(nodesCoord(:,1),xi);
            d = dxydeta(nodesCoord(:,2),xi);
            nodalJac(i,k) = (a*d - b*c)^2;
            k = k + 1;
        end
    end
end

elJacRatio = min(nodalJac,[],2)./max(nodalJac,[],2);

minJacRatioMesh = min(elJacRatio);

if minJacRatioMesh > 0.5
    fprintf("\nMinimum Jacobian ratio in the mesh = %.2f\n",minJacRatioMesh)
else
    warning("[warning] Minimum Jacobian ratio in the mesh is < 0.5, highly distorted" + ...
        " elements!")
    numEldistorted = sum(elJacRatio < 0.5);
    percDistorted = numEldistorted/size(elements,1)*100;
    fprintf("\nNumber of highly distorted elements: %d (%.2f%%)\n",numEldistorted, ...
        percDistorted)
end

if plotBool
    figure
    histogram(elJacRatio,"Normalization","percentage")
    grid on
    ylabel("Frequency [%]","FontSize",18)
    xlabel("Jacobian ratio (best close to 1)","FontSize",18)
    title("Mesh Jacobian ratio distribution","FontSize",21)
end

end

