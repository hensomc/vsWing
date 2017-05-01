function epsRitz=getRitzStrain( solDATA, panXY, lamDATA, surfI, C, W,dWx,dWy,dWxy,dWxx,dWyy, psiMesh, etaMesh )

% Recovers the Ritz strain components given the derivative functions Fx,
% Fy, Fxy and the Ritz coeff vector C.

% Strains are recovered on a grid of points: [psiMesh, etaMesh] which must be
% consistent with the evaluations of Fx, Fy, Fxy.

% Recover solution parameters
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

zvec=zeros(1,(pDeg+1)^2);
dim=(pDeg+1)^2;

% Matrix             dimensions
% =================  ===========================
% psiMesh, etaMesh:  [(Mg+1)^2,     1]
% U, V, W, P, Q, R:  [(Mg+1)^2,    5*(nDeg+1)^2]
%                C:  [5*(nDeg+1)^2, 1];

% Displacement basis functions:
U=W;dUx=dWx;dUy=dWy;dUxy=dWxy;dUxx=dWxx;dUyy=dWyy;
V=W;dVx=dWx;dVy=dWy;dVxy=dWxy;dVxx=dWxx;dVyy=dWyy;
P=W;dPx=dWx;dPy=dWy;dPxy=dWxy;dPxx=dWxx;dPyy=dWyy;
R=W;dRx=dWx;dRy=dWy;dRxy=dWxy;dRxx=dWxx;dRyy=dWyy;

% outer fiber z-position for strain calcs @ outer fiber
z=sum(lamDATA.thk)/2;


% Reference mid-surface
surfP=[];
for i=1:size(surfI,1)
    surfP(i).xyz=getWingZCoords(surfI(i),psiMesh,etaMesh,panXY,i-1,0);  % integ pts on wing surface
end


epsRitz=[];
% Calc strain tensor at mesh points: length=(Mg+1)^2
for i=1:length(psiMesh)

    % Offset from reference surface to laminate midplane
    z0=0;
    if ~isempty(surfP)
        z0=surfP(1).xyz(i,3);
    end
    
    %Strain basis functions: [1, 5*(nDeg+1)^2]
    Ux=[dUx(i,:) zvec zvec zvec zvec]*C;
    Uy=[dUy(i,:) zvec zvec zvec zvec]*C;
    Vx=[zvec dVx(i,:) zvec zvec zvec]*C;
    Vy=[zvec dVy(i,:) zvec zvec zvec]*C;
    Wx=[zvec zvec dWx(i,:) zvec zvec]*C;
    Wy=[zvec zvec dWy(i,:) zvec zvec]*C;
    Px=[zvec zvec zvec dPx(i,:) zvec]*C;
    Py=[zvec zvec zvec dPy(i,:) zvec]*C;
    Rx=[zvec zvec zvec zvec dRx(i,:)]*C;
    Ry=[zvec zvec zvec zvec dRy(i,:)]*C;
    r =[zvec zvec zvec P(i,:)   zvec]*C;
    p =[zvec zvec zvec zvec   R(i,:)]*C;
    
    % Strain matrix [epsbar]: [12, 5*(nDeg+1)^2]
    epsbar=[Ux; Uy; Vx; Vy; Wx; Wy; Px; Py; Rx; Ry; p; r];

    % Jacobian
    J=jacob2D_Iso( psiMesh(i), etaMesh(i), panXY );
    JI=inv(J);
    j11=JI(1,1); j12=JI(1,2); j21=JI(2,1); j22=JI(2,2);
    
    % Strain transformation matrix: [5, 12]
    T0=[j11 j12  0   0   0   0   0    0    0    0    0   0;
        0   0  j21 j22  0   0   0    0    0    0    0   0;
        j21 j22 j11 j12  0   0   0    0    0    0    0   0;
        0   0   0   0  j21 j22  0    0    0    0    0   1;
        0   0   0   0  j11 j12  0    0    0    0    1   0];
    
    T1=[ 0   0   0   0   0   0  j11  j12   0    0    0   0;
        0   0   0   0   0   0   0    0   j21  j22   0   0;
        0   0   0   0   0   0  j21  j22  j11  j12   0   0;
        0   0   0   0   0   0   0    0    0    0    0   0;
        0   0   0   0   0   0   0    0    0    0    0   0];

    if(lamDATA.tDeg > 0)
        z=0.5*thk_legendre(lamDATA.tDeg,psiMesh(i), etaMesh(i))*lamDATA.thk_coeff'; % assume symmetric distribution
    end
    
    % Compute strain tensor: [5 , 12]*[12, 5*(nDeg+1)^2] = [5,1]
    epsRitz=[epsRitz ( T0 + (z0+z)*T1)*epsbar];  % epsRitz = [epsxx epsyy epsxy epsyz epsxz]
end

% Transpose to get [(Mg+1)^2, 5];
epsRitz=epsRitz';

