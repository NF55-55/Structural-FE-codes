clear
close all
clc

tic
% profile on
%% Quick simulation controls. For the fondamental ones, like mat. prop. scroll down
plotMesh            = false;% to visualize nodes numbers, available only for small mesh
checkMesh           = false;% checking mesh quality
plotJacRatioDistr   = false;% this is plotted only in case of highgly distorted elem.
t                   = 2;    % mm (part/element thickness)
planeStress         = true; % if false then plane strain is assumed with infinite t
visualizeMsparsity  = false;% visualize the sparsity of [M_global] matrix
dt                  = 1e-2; % s initial time step size
t0                  = 0;    % s start time
tFinal              = 1;    % s duration of simulation
gamma               = 0.5;  % constant of the Newmark-beta implicit time integration
beta                = 0.25; % constant of the Newmark-beta implicit time integration
loadFromStatic      = false;% load initial displacement field from static analysis
maxNewtonIters      = 5;    % maximum number of Newton-Raphson loops per time step
tolNewton           = 5e-3; % tolerance used for convergence of Newton method (residual)
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
         0      0       1/G];

    S = inv(C); % stiffness matrix - material constitutive law
else
    C = ((1+nu)/E)*[(1-nu)   -nu   0;
                    -nu    (1-nu)  0;
                    0      0       2];
    S = inv(C);
end
%% Indices of Boundary Conditions and Loads

% d.o.f. or nodal displacement vector [mm]
d = nan(2*size(nodes,1),1); % {d11; d12; d21; d22; d31; d32; d41; ...}
% where first subscript is the node and second one is the direction, they
% are initialized to NaN so that is possible to distinguish the unknown
% d.o.f. from the one known from BCs.

% nodal force vector [N]
Fref = zeros(size(nodes,1)*2,1);

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

Fref(~isnan(d)) = nan; % only the constrained dof set to nan

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

%% Solve implicit NON-linear dynamics (Newmark-beta)

% Initial conditions I.C.
% Creating the problem variables for old time instant t_k = t0(in this case)
if loadFromStatic
    load("..\..\static2dynamics\dispField.mat","dFromStatic")
    fprintf("\nDisplacement field imported from static solver as requested.\n")
    dk = dFromStatic;
else
    dk = zeros(size(d));
end
velk = zeros(size(d)); % nodal velocity vector [mm/s]
Fext = zeros(size(d)); % also initial condition
Fint = FintCompute(S,nodes,elements,dk,t);
acck = M \ (Fext - Fint); % nodal acceleration vector [mm/s^2]

% Amplitude of load
Fampx = -3000;
Fampy = -Fampx;

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

fprintf("\nSolving NON-LINEAR implicit dynamics...\n")

while (time < tFinal+dt)
    % Imposing forces at time instant t_{k+1}
    Fext(dof2index(FxNodes,'x')) = Fampx/(length(FxNodes))*time/tFinal;
    Fext(dof2index(FyNodes,'y')) = Fampy*time/tFinal;
    Fext = rearrangeFext(Fref,Fext);
    % initializing the nodal displacement vector for new time instant
    % dk1 = dk                                  -> more stable but maybe slower
    % dk1 = dk + velk*dt + (0.5-beta)*acck*dt^2 -> less stable but faster
    dk1 = dk;
    % Imposing boundary conditions at time instant t_{k+1}
    dk1(dof2index(BCxNodes,'x')) = 0;
    dk1(dof2index(BCyNodes,'y')) = 0;

    % Initializing the maximum normalized residual force
    normR  = 100;
    % Newton-Raphson iterations counter
    newtonCounter = 0;

    % Newton-Raphson loop
    while (normR > tolNewton)
        % computing F^{int} based on the current u^{k+1}
        % as well as computing KtMat and KtGeo (for efficiency)
        [Fint, KtMat, KtGeo] = FintAndKtInt(S,nodes,elements,dk1,t);
        % no external force stiffness matrix as we assume load direction
        % and magnitude is independent from the displacement field
        % Effective tangent stiffness matrix for the current u^{k+1}:
        KtEff = (1/(beta*dt^2))*M + KtMat + KtGeo;
        % computing nodal acceleration vector at t^{k+1} (Newmark-beta)
        acck1 = (1/(beta*dt^2))*(dk1 - dk - velk*dt - (0.5-beta)*acck*dt^2);
        % computing nodal velocity vector at t^{k+1} (Newmark-beta)
        velk1 = velk + (1-gamma)*acck*dt + gamma*dt*acck1;
        % Split M, Fint and Fext into known forces and unknown ones
        % so to exclude from computations the unknown reaction forces
        [M,Fint,acck1,ordRow,ordCol,ind] = ...
                            rearrangeUnknowns(M,Fref,d,Fint,acck1);
        % ind is used to identify the first ind terms of Fint and Fext
        % which are not related to reaction forces (unknowns)
        % computing current (t^{k+1}) equilibrium residual
        R = M(1:ind,1:ind)*acck1(1:ind) + Fint(1:ind) - Fext(1:ind);
        % For the convergence criterion
        normR = max(abs(R)/max(abs(Fext)));%norm(R);
        % Rearranging columns so to only compute the unknown delta(d) for
        % the unknown dof, the rest of them are assumed to have delta(d) =
        % 0 as they're constrained
        [KtEff,dk1,~] = rearrangeUnknowns2(Fref,d,KtEff,dk1);
        
        if (~warnCondNum)
            % check conditioning number before inversion of Kt effective:
            condKtEff = condest(KtEff(1:ind,1:ind));
            if (condKtEff > maxConditionNumb)
                warning("[Warning] The matrix to be inverted is not well conditioned!")
                fprintf("\nCondition number of KtEff(1:ind,1:ind) = %.2e\n",condKtEff)
            end
            warnCondNum = true;
        end
        % solving the linearized problem R = R_{k+1} + dR/du|_{k+1}*delta(u) = 0
        % finding delta(u) = - (dR/du|_{k+1})^{-1} * R_{k+1}
        % so the new u_{k+1} = u_{k} + delta(u):
        dk1(1:ind) = dk1(1:ind) - KtEff(1:ind,1:ind) \ R;
        [~, M] = ordRows(ordRow,M);
        [~, M] = ordCols(ordCol,M);
        [~, dk1] = ordCols(ordCol,dk1);
        % Incrementing the Newton-Raphson iter counter
        newtonCounter = newtonCounter + 1;
        % stop simulation if number of newton iters is larger than the
        % amount set by user
        if (newtonCounter > maxNewtonIters)
            error("[Error] Maximum number of Newton-Raphson loops reached!")
        end
    end
    [~, acck1] = ordCols(ordCol,acck1);
    [~, Fext] = ordRows(ordRow,Fext);
    fprintf("Newton-Raphson iterations: %d\n",newtonCounter)
    
    % Saving results in vtk format:
    if (time >= tSaveVTK)
        vtkResultsGen(nodes, elements, dk1, S, countVTK, time)
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