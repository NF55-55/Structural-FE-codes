function matlabPlotResults(nodes, elements, d, S, defScaleFact)

% plot displacement field (linear node interpolation as in commercial FEM
% codes):
dplot = d*defScaleFact;
dMag = sqrt(d(1:2:end-1).^2 + d(2:2:end).^2);

% contour plot of x-displacement field
figure
colormap("jet")
for i = 1:size(elements,1)
    patch(nodes(elements(i,:),1) + dplot(dof2index(elements(i,:),'x')), ...
        nodes(elements(i,:),2) + dplot(dof2index(elements(i,:),'y')), ...
        d(dof2index(elements(i,:),'x')),'FaceColor', 'interp', ...
        'EdgeColor', 'black')
    hold on
end
title("Contour plot of d_x field with deformation scale factor = "+ ...
    num2str(defScaleFact))
axis equal
grid on
colorbar;
hold off

% contour plot of y-displacement field
figure
colormap("jet")
for i = 1:size(elements,1)
    patch(nodes(elements(i,:),1) + dplot(dof2index(elements(i,:),'x')), ...
        nodes(elements(i,:),2) + dplot(dof2index(elements(i,:),'y')), ...
        d(dof2index(elements(i,:),'y')),'FaceColor', 'interp', ...
        'EdgeColor', 'black')
    hold on
end
title("Contour plot of d_y field with deformation scale factor = "+ ...
    num2str(defScaleFact))
axis equal
grid on
colorbar;
hold off

% contour plot of displacement magnitude field
figure
colormap("jet")
for i = 1:size(elements,1)
    patch(nodes(elements(i,:),1) + dplot(dof2index(elements(i,:),'x')), ...
        nodes(elements(i,:),2) + dplot(dof2index(elements(i,:),'y')), ...
        dMag(elements(i,:)),'FaceColor', 'interp', ...
        'EdgeColor', 'black')
    hold on
end
title("Contour plot of d_{Mag} field with deformation scale factor = "+ ...
    num2str(defScaleFact))
axis equal
grid on
colorbar;
hold off

% Stresses visualization

figure
colormap("jet")
for i = 1:size(elements,1)
    dTemp = d(2*elements(i,:)+[-1;0]);
    sigmax = elStresses(nodes(elements(i,:),:),S,reshape(dTemp,[],1),'x');
    sigmax = extrapSigmaNodes(sigmax);
    patch(nodes(elements(i,:),1) + dplot(dof2index(elements(i,:),'x')), ...
        nodes(elements(i,:),2) + dplot(dof2index(elements(i,:),'y')), ...
        sigmax,'FaceColor', 'interp', ...
        'EdgeColor', 'black')
    hold on
end
title("Contour plot of sigma_x field with deformation scale factor = "+ ...
    num2str(defScaleFact))
axis equal
grid on
colorbar;
hold off 

end

