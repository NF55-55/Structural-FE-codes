function [Fint, KtMat, KtGeo] = FintAndKtInt(S,nodes,elements,dk1,t)

% As Fint, KtMat and KtGeo all require Gauss quadrature over the elements
% and use F and B0, it's much more efficient to compute them together in
% one single loop over the elements instead of using 3 of them
% Inputs:   - S: material constitutive stiffness matrix
%           - nodes: coordinates of each node
%           - elements: connectivity table (nodeIDs per element)
%           - dk1: current nodal displacement vector d_{k+1}
%           - t: elements thickness

% problem size
sizProb = length(dk1);
nEl     = size(elements,1);

Fint = zeros(sizProb,1); % initializing the output vector

% same for both KtGeo and KtMat
ivect = zeros(nEl*64,1);
jvect = zeros(nEl*64,1);

valVectMat = zeros(nEl*64,1);

valVectGeo = zeros(nEl*64,1);

w = 1; % gaussian weights

for i = 1:nEl
    % indices of element's i dof
    dofID = element2Ind(elements(i,:));
    % current nodal disp. vector for element i
    dk1El = dk1(dofID);
    % Lagrangian coordinates
    nodesCoord0 = nodes(elements(i,:),:)';
    % current nodal coordinates
    eulerCoord = nodesCoord0(:) + dk1El;
    % re-initializing KtMatEl
    KtMatEl = zeros(numel(eulerCoord)*[1,1]);
    % re-initializing KtGeoEl
    KtGeoEl = zeros(length(dk1El),length(dk1El));
    % looping through the gauss points of the i-th element
    for xi=[-1/sqrt(3), 1/sqrt(3)]
        for eta=[-1/sqrt(3), 1/sqrt(3)]
            % material derivatives of the shape functions [B0] and
            % displacement gradient tensor [F] in material coordinates X_i
            [B0,F,detJac] = FandB0computeOptim(nodesCoord0',xi,eta,eulerCoord);
            % Green-Lagrange strain tensor
            E = 1/2 * (F' * F - eye(size(F)));
            % Green-Lagrange strain tensor (Voigt notation)
            E = [E(1,1); E(2,2); 2*E(1,2)];
            % Piola-Kirchhoff second stress tensor (Voigt notation)
            PK2 = S * E;
            % Piola-Kirchhoff second stress tensor
            PK2 = [PK2(1), PK2(3);
                   PK2(3), PK2(2)];
            % Nominal stress tensor
            P   = F * PK2;
            % Nominal stress in vectorial format: [Pxx; Pxy; Pyx; Pyy] see
            % notes (to be consistent with the derivative matrix)
            P = [P(1,1); P(1,2); P(2,1); P(2,2)];
            % Fint computation---------------------------------------------
            Fint(dofID) = Fint(dofID) + (t * w*w * B0' * P * detJac);
            % KtMatEl computation------------------------------------------
            % computing [Edot*] which is in Edot = [Edot*]udot
            % it should be in the order [E11;E12;E22], but temporarely it will
            % be [E11,E22,E12] to be compatible with [C] (here called S)
            EdotStar = [(B0(1,:)*F(1,1) + B0(3,:)*F(2,1));
                        (B0(2,:)*F(1,2) + B0(4,:)*F(2,2));
                    (B0(1,:)*F(1,2) + B0(3,:)*F(2,2) + B0(2,:)*F(1,1) + B0(4,:)*F(2,1))];
            SdotStar = S * EdotStar; % still in the "wrong" order [S11,S22,S12]
            % re-order (and add) SdotStar rows so to give the order compatible
            % with B0, F, ... which is [S11,S12,S21,S22]
            SdotStar = [SdotStar(1,:);
                        SdotStar(3,:);
                        SdotStar(3,:); % S21 = S12
                        SdotStar(2,:)];
            % now we compute PdotStar
            PdotStar = [(F(1,1)*SdotStar(1,:) + F(1,2)*SdotStar(3,:));
                        (F(1,1)*SdotStar(2,:) + F(1,2)*SdotStar(4,:));
                        (F(2,1)*SdotStar(1,:) + F(2,2)*SdotStar(3,:));
                        (F(2,1)*SdotStar(2,:) + F(2,2)*SdotStar(4,:))];
            KtMatEl = KtMatEl + t*w*w*(B0' * PdotStar)*detJac;
            % KtGeoEl computation------------------------------------------
            % computing [Fdot][S]* that is [Fdot][S] = ([Fdot][S]*)udot
            FdotSstar = [(B0(1,:)*PK2(1,1) + B0(2,:)*PK2(2,1));
                         (B0(1,:)*PK2(1,2) + B0(2,:)*PK2(2,2));
                         (B0(3,:)*PK2(1,1) + B0(4,:)*PK2(2,1));
                         (B0(3,:)*PK2(1,2) + B0(4,:)*PK2(2,2))];
            % KtGeo for element i at the (xi,eta) integration/Gauss point
            KtGeoEl = KtGeoEl + t*w*w*(B0' * FdotSstar)*detJac;
        end
    end
    % for sparse matrix creation
    iTemp = dofID*ones(1,8);
    ivect(i*64-63:i*64) = iTemp(:);
    jTemp = ones(8,1)*dofID';
    jvect(i*64-63:i*64) = jTemp(:);
    % For KtMatEl
    valVectMat(i*64-63:i*64) = KtMatEl(:);
    % For KtGeoEl
    valVectGeo(i*64-63:i*64) = KtGeoEl(:);
end

% Assemblying the sparse matrices

% KtMat
KtMat = sparse(ivect,jvect,valVectMat,sizProb,sizProb);
% KtGeo
KtGeo = sparse(ivect,jvect,valVectGeo,sizProb,sizProb);

end

