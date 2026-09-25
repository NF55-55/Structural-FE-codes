clear
close all
clc
% profile on
tic
%% Quick simulation controls. For the fondamental ones, like mat. prop. scroll down
plotMesh            = false; % to visualize nodes numbers, available only for small mesh
checkMesh           = false; % checking mesh quality
plotJacRatioDistr   = false; % this is plotted only in case of highgly distorted elem.
t                   = 2;     % mm (part/element thickness)
planeStress         = true;  % if false then plane strain is assumed with infinite t
visualizeKsparsity  = false; % visualize the sparsity of [K_global] matrix
%defScaleFact        = 5;     % scale factor used in the plots
%typeStress2Plot     = 'x';   % available are 'x', 'y', 'xy' (tau), 'all'
transferTOdynamics  = false;  % transfer displacement field to initialize dynamic sim.

%% Import mesh nodes and elements from .inp file of Abaqus or Ansys
% nodes of each element should be already in CCW numbering
[nodes, elements] = importMesh();

% plot the geometry and mesh
if plotMesh && size(nodes,1)<5000
    figure
    for i = 1:size(elements,1)
        plot(nodes([elements(i,:) elements(i,1)],1),nodes([elements(i,:) ...
            elements(i,1)],2),'k-*')
        hold on
    end
    grid on
    axis equal
    title("2D Mesh","FontSize",18)
    xlabel("x direciton [mm]","FontSize",15)
    ylabel("y direction [mm]","FontSize",15)
    for i = 1:size(nodes,1)
        text(nodes(i,1),nodes(i,2),num2str(i),"FontSize",12,"Color",[1,0,0.1], ...
            "VerticalAlignment","bottom","HorizontalAlignment","left")
    end
    hold off
    tempVar = input("Take your time writing the BCs and loads looking at" + ...
        " the nodes numbers, press [any number]+[enter] to terminate.");
    return
end

%% Check mesh quality (Jacobian ratio)
if checkMesh
    checkMeshJac(nodes,elements,plotJacRatioDistr)
end

%% Material properties (homogeneous-isotropic-elastic)

E = 200e3; % MPa (Young's modulus)
nu = 0.33; % (Poisson's ratio)
G = E/(2*(1+nu)); % MPa (Shear modulus)

if planeStress
    C = [1/E    -nu/E   0;
         -nu/E  1/E     0;
         0      0       1/G]; % compliance matrix - material constitutive law

    S = inv(C); % stiffness matrix - material constitutive law
else
    S = E/((1+nu)*(1-2*nu))*[1-nu,nu,0; nu,1-nu,0; 0,0,(1-2*nu)/2];
end
%% Boundary Conditions and Loads

% d.o.f. or nodal displacement vector [mm]
d = nan(2*size(nodes,1),1); % {d11; d12; d21; d22; d31; d32; d41; ...}
% where first subscript is the node and second one is the direction, they
% are initialized to NaN so that is possible to distinguish the unknown
% d.o.f. from the one known from BCs.

% nodal force vector [N]
F = zeros(size(nodes,1)*2,1);

% apply boundary conditions (BCs):
BCxNodes = load(".\BC and load node sets\BCx.txt");
BCyNodes = load(".\BC and load node sets\BCy.txt");
FxNodes = load(".\BC and load node sets\Fx.txt");
FyNodes = load(".\BC and load node sets\Fy.txt");

if ~isempty(BCxNodes)
    d(dof2index(BCxNodes,'x')) = 0;
end
if ~isempty(BCyNodes)
    d(dof2index(BCyNodes,'y')) = 0;
end

% apply loads
if ~isempty(FxNodes)
    F(dof2index(FxNodes,'x')) = 10000/length(FxNodes);
end
if ~isempty(FyNodes)
    F(dof2index(FyNodes,'y')) = 0;
end

% rearrange [K], {d} and {F} so to have at the top the unknown
% displacements and at the bottom the equations referring to the unknown
% forces. Also store the new order so to be able to restore {d} and {F} at 
% the end of computations to original order to plot the results.
[d,F,mapVect] = rearrangeUnknowns(d,F);

%% Global Stiffness Matrix calculation

%profile on

% Matlab's sparse(i,j,v,size1,size2) function to create a sparse matrix automatically
% adds up the entries in the same element of the matrix. i and j are the
% rows and cols for the values v to enter, while size1 and size2 are the
% size specifiers for the global K matrix. To be more efficient all the
% elements of K_el are stored in the following vectors together with the
% "coordinates" in the K matrix. After having computed K_el of every
% element, the K sparse matrix will be created.

ivect = zeros(size(elements,1)*64,1);
jvect = zeros(size(elements,1)*64,1);
valVect = zeros(size(elements,1)*64,1);

nel1perc = size(elements,1)/100;
nMessage = nel1perc;

for i=1:size(elements,1)
    nodalCoords = nodes(elements(i,:),:);
    Kel = elStiffMat(nodalCoords,S,t);
    dofNumb = [2*elements(i,:)-1; 2*elements(i,:)];
    dofNumb = dofNumb(:);
    valVect(i*64-63:i*64) = Kel(:);
    iTemp = dofNumb*ones(1,8);
    ivect(i*64-63:i*64) = iTemp(:);
    jTemp = ones(8,1)*dofNumb';
    jvect(i*64-63:i*64) = jTemp(:);
    if i>nMessage
        fprintf("\nComputing global stiffness matrix... %.2f%% done.", ...
            nMessage/nel1perc)
        nMessage = nMessage + nel1perc;
    end
end
K = sparse(ivect,jvect,valVect,length(d),length(d));
K = K(mapVect, mapVect); % setting K in the solving order
clear ivect jvect valVect
fprintf("\nGlobal stiffness matrix [K] computed successfully!\n")

% To understand the huge computational cost of full matrix K, for a small mesh:
% Time required with full matrix: 316.08s. Memory usage for [K]: 238.5   MB
% Time required with sparse matrix: 3.79s. Memory usage for [K]:   0.045 MB

% plot [K] sparsity
if visualizeKsparsity
    figure
    spy(K)
    title("Global Stiffness Matrix [K] sparsity visualization","FontSize",21)
end

%% Compute the static equilibrium - FE main output {d}

fprintf("\nSolving the system of linear algebraic equations [K]{u}={F}...\n")

% solve for unknown displacements
% First, find the n° of rows to consider
ind = sum(isnan(d)); % length of {d_unkn}

% solve {d_unkn}:

% now the system is: [k1 , k2]{u(unknown); u(known)} = {F(known)}
% so: [k1]*{u(unk)} + [k2]*{u(kno)} = {F(kno)}
d(1:ind) = K(1:ind,1:ind)\(F(1:ind)-K(1:ind,ind+1:end)*d(ind+1:end));

fprintf("\nUnknown displacement computed successfully!\n")

%F(ind+1:end) = K(ind+1:end,:)*d;
%fprintf("\nReaction forces computed successfully!\n")

d(mapVect) = d; % {d} final (in the original order)

%F(mapVect) = F; % {F} final (in the original order)

% profile off
% profile viewer

%% Displacements visualization and transfer displacement field for dynamics

% for large meshes it's not feasible in MATLAB, logical statement required to check
% problem size automatically, if too large, create .vtk file, so to
% visualize results in Paraview

% if size(nodes,1)<100
%     matlabPlotResults(nodes, elements, d, S, defScaleFact)
% else

% Actually, it's always better to save the results in .vtk and then visualize
% and do post-processing in Paraview
vtkResultsGen(nodes, elements, d, S)

if transferTOdynamics
    dFromStatic = d;
    save("..\static2dynamics\dispField","dFromStatic")
    fprintf("\nDisplacement field made available to the dynamic solver as requested.\n")
end

toc
% profile off
% profile viewer