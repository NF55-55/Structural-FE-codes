clear
close all
clc

% profile on
tic
%% Quick simulation controls. For the fondamental ones, like mat. prop. scroll down
plotMesh            = false;% to visualize nodes numbers, available only for small mesh
checkMesh           = false;% checking mesh quality
plotJacRatioDistr   = false;% this is plotted only in case of highgly distorted elem.
t                   = 2;    % mm (part/element thickness)
planeStress         = true; % if false then plane strain is assumed with infinite t
visualizeMsparsity  = false;% visualize the sparsity of [M_global] matrix
visualizeKsparsity  = false;% visualize the sparsity of [K_global] matrix
dt                  = 1e-2; % s initial time step size
t0                  = 0;    % s start time
tFinal              = 1;    % s duration of simulation
gamma               = 0.5;  % constant of the Newmark-beta implicit time integration
beta                = 0.25; % constant of the Newmark-beta implicit time integration
loadFromStatic      = false; % load initial displacement field from static analysis
timeInstants2save   = 50;   % uniformly distributed time instants to save as .vtk files
maxConditionNumb    = 1e10; % maximum condition number of the effective tang stiff mat
                            % knowing that lost significant digits are proportional to
                            % the order of mag. of condition number. For double
                            % precision significant digits are 17 - LOG10(ConditiNumb);
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
    return
end

%% Check mesh quality (Jacobian ratio)
if checkMesh
    checkMeshJac(nodes,elements,plotJacRatioDistr)
end

%% Material properties (homogeneous-isotropic-elastic)

rho = 7.85e-9; % t/mm^3

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
%% Indices of Boundary Conditions and Loads

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

dRef = d;

% apply loads (NO LOADS APPLIED AT t0)
% if ~isempty(FxNodes)
%     F(dof2index(FxNodes,'x')) = 2e3/length(FxNodes);
% end
% if ~isempty(FyNodes)
%     F(dof2index(FyNodes,'y')) = 0;
% end

F(~isnan(d)) = nan; % only the constrained dof set to nan

%% Global Mass Matrix calculation

% Matlab's sparse(i,j,v,size1,size2) function to create a sparse matrix automatically
% adds up the entries in the same element of the matrix. i and j are the
% rows and cols for the values v to enter, while size1 and size2 are the
% size specifiers for the global M matrix. To be more efficient all the
% elements of M_el are stored in the following vectors together with the
% "coordinates" in the M matrix. After having computed M_el of every
% element, the M sparse matrix will be created.
ivect = zeros(size(elements,1)*64,1);
jvect = zeros(size(elements,1)*64,1);
valVect = zeros(size(elements,1)*64,1);

nel1perc = round(size(elements,1)/100);
nMessage = nel1perc;

for i=1:size(elements,1)
    nodalCoords = nodes(elements(i,:),:);
    Mel = elMassMat(nodalCoords,rho,t);
    dofNumb = [2*elements(i,:)-1; 2*elements(i,:)];
    dofNumb = dofNumb(:);
    valVect(i*64-63:i*64) = Mel(:);
    iTemp = dofNumb*ones(1,8);
    ivect(i*64-63:i*64) = iTemp(:);
    jTemp = ones(8,1)*dofNumb';
    jvect(i*64-63:i*64) = jTemp(:);
    if i>nMessage
        fprintf("\nComputing global Mass matrix... %.2f%% done.", ...
            i/size(elements,1)*100)
        nMessage = nMessage + nel1perc;
    end
end
M = sparse(ivect,jvect,valVect,length(d),length(d));
clear ivect jvect valVect
fprintf("\nGlobal Mass matrix [M] computed successfully!\n")

% plot [M] sparsity
if visualizeMsparsity
    figure
    spy(M)
    title("Global Stiffness Matrix [K] sparsity visualization","FontSize",21)
end

%% Global Stiffness Matrix calculation

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

nel1perc = round(size(elements,1)/100);
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
            i/size(elements,1)*100)
        nMessage = nMessage + nel1perc;
    end
end
K = sparse(ivect,jvect,valVect,length(d),length(d));
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

%% Solve implicit linear dynamics (Newmark-beta)

% rearrange [K], {d} and {F} so to have at the top the unknown
% displacements and at the bottom the equations referring to the unknown
% forces. Also store the order of the rows and cols of K so to be able to
% restore {d} and {F} at the end of computations to original order, to be
% able to plot the results.
[M,K,d,ordRow,ordCol] = rearrangeUnknowns(M,K,d,F);

% First, find the n° of rows to consider
ind = sum(isnan(d)); % length of {d_unkn}

% The system of equations is described in "Newmark-beta and Newton-Raphson"
% notes, first we solve for the unknown displ. and then for the reaction
% forces, just like in statics, but now using Newmark-beta implicit time
% integration scheme.

% Initial conditions I.C.
% Creating the problem variables for old time instant t_k = t0(in this case)
if loadFromStatic
    load("..\..\static2dynamics\dispField.mat","dFromStatic")
    fprintf("\nDisplacement field imported from static solver as requested.\n")
    % we need to re-arrange dFromStatic into the current order (unk-kno)
    dTemp = dFromStatic;
    dTemp = [dTemp(isnan(dRef)); dTemp(~isnan(dRef))];
    dk = dTemp;
else
    dk = zeros(size(d));
end
velk = zeros(size(d)); % nodal velocity vector [mm/s]
% Initializing the problem variables for new time instant t_{k+1}
dk1   = zeros(size(d));
velk1 = zeros(size(d));
acck1 = zeros(size(d));
Fk1   = zeros(size(d)); % also initial condition
acck = M \ (Fk1 - K*dk); % nodal acceleration vector [mm/s^2]

% Amplitude of load
Fampx = 3000;
Fampy = 3000;

% Time step for saving:
dtSaveVTK = tFinal/timeInstants2save;
tSaveVTK  = t0 + dtSaveVTK;
countVTK  = 1;

% Time variable
time = t0 + dt;

% deletes all .vtk files (if there's any) in the VTKresults folder
dirTemp = dir(".\VTK RESULTS\");
if contains([dirTemp.name],'.vtk')
    fprintf("\nDeleting the results in ./VTK RESULTS folder...\n")
    delete(".\VTK RESULTS\*.vtk");
end
clear dirTemp;

% Warning about condition number of KtEff initialization
warnCondNum = false;

fprintf("\nSolving implicit dynamics...\n")

while (time < tFinal+dt)
    % Imposing forces at time instant t_{k+1}
    % Ftemp is initially referred to as in the original order of the d.o.f.
    Ftemp = zeros(size(Fk1));
    % assigning the forces to the specific dof chosen with the sets
    Ftemp(dof2index(FxNodes,'x')) = Fampx/(length(FxNodes))*time/tFinal;
    Ftemp(dof2index(FyNodes,'y')) = Fampy*time/tFinal;
    % ordering this Ftemp so to match the Fk1 order
    Ftemp = [Ftemp(~isnan(F)); Ftemp(isnan(F))];
    % assignment of forces to Fk1
    Fk1 = Ftemp;
    % Imposing boundary conditions at time instant t_{k+1}
    % same as in Fk1:
    dxTemp = nan(size(d));
    dxTemp(dof2index(BCxNodes,'x')) = 0;
    dxTemp(dof2index(BCyNodes,'y')) = 0;
    dxTemp = [dxTemp(isnan(dRef)); dxTemp(~isnan(dRef))];
    dk1(~isnan(dxTemp)) = dxTemp(~isnan(dxTemp));
    
    % Solving for the unknown nodal displacements d_{k+1}^{unk}
    % check conditioning number before inversion:
    if (~warnCondNum)
        % check conditioning number before inversion of K effective:
        condKEff = condest(1/(beta*dt^2)*M(1:ind,1:ind) + K(1:ind,1:ind));
        if (condKEff > maxConditionNumb)
            warning("[Warning] The matrix to be inverted is not well conditioned!")
            fprintf("\nCondition number of KtEff(1:ind,1:ind) = %.2e\n",condKtEff)
        end
        warnCondNum = true;
    end
    dk1(1:ind) = (1/(beta*dt^2)*M(1:ind,1:ind) + K(1:ind,1:ind)) \ (1/(beta*dt^2) ...
        *M(1:ind,1:ind)*(dk(1:ind) + velk(1:ind)*dt + (0.5-beta)*acck(1:ind)*dt^2) + ...
        1/(beta*dt^2)*M(1:ind,ind+1:end)*(-dk1(ind+1:end) + dk(ind+1:end) + ...
        velk(ind+1:end)*dt + (0.5-beta)*acck(ind+1:end)*dt^2) - K(1:ind,ind+1:end) ...
        *dk1(ind+1:end) + Fk1(1:ind));
    % Solve for velocity vel_{k+1}
    velk1 = (gamma/(beta*dt))*(dk1 - dk) + (1-gamma/beta)*velk + (1-gamma/(2*beta)) ...
        *dt*acck;
    acck1 = 1/(beta*dt^2)*(dk1 - dk - velk*dt - (0.5-beta)*acck*dt^2);
    % computation of reaction forces is not necessary if not explicitly requested
    % if so, just solve the lower part of the system of eqs (rows: ind+1:end)
    
    % Saving results in vtk format:
    if (time >= tSaveVTK)
        % providing displacements in the original order
        [~, dk1Rearr] = ordCols(ordCol,dk1);
        vtkResultsGen(nodes, elements, dk1Rearr, S, countVTK, time)
        tSaveVTK = tSaveVTK + dtSaveVTK;
        countVTK = countVTK + 1;
        fprintf("\nTime instant %f / %f s (%.2f %%) saved.\n",time,tFinal, ...
            time/tFinal*100)
    end

    % updating the variables vector for next time instant
    dk   = dk1;
    velk = velk1;
    acck = acck1;

    time = time + dt;
end

fprintf("\nSimulation completed!\n")

fprintf("\nTo display time in paraview: filter -> python annotation -> " + ...
    "'Time: %%.3e' %% (input.FieldData['TIME'][0])\n")

elapsedTime = toc;
fprintf("\nSimulation runtime: %d minutes and %.1f s\n",floor(elapsedTime/60),...
    elapsedTime-floor(elapsedTime/60)*60)
% profile off
% profile viewer