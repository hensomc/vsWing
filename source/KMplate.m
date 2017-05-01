function [K,M,Wt]=KMplate(pDeg,intMethod,smearKM,panXY,lamDATA,psiI,etaI,zI,gl_wts,B,plot_flag)

if(plot_flag>0)
    disp('Calculating stiffness matrix...');
end
ndof=5*(pDeg+1)^2;
K=zeros(ndof,ndof);M=zeros(ndof,ndof);
zvec=zeros(1,(pDeg+1)^2);

Wt=0;
A=0;

nply = length(lamDATA.thk);

% Layer thickness functions
thk=lamDATA.thk;

psi_current = -9999;

% Recover basis vectors from [B]
U=B.U;
dUx=B.Ux;
dUy=B.Uy;

V=B.V;
dVx=B.Vx;
dVy=B.Vy;

W=B.W;
dWx=B.Wx;
dWy=B.Wy;

P=B.P;
dPx=B.Px;
dPy=B.Py;

R=B.R;
dRx=B.Rx;
dRy=B.Ry;

% Transverse shear correction coefficient
%Kshear = 5/6;  % For thick plates Kshear = 5/6 or Pi^2/12 and 
Kshear = 0.025;

for i=1:length(psiI);
    gl_wt=1.0;
    if(intMethod == 2)  %gaussian quadrature weights (in 2D)
        gl_wt=gl_wts(i);
    end
    
    % Strain basis functions
    Ux=[dUx(i,:) zvec zvec zvec zvec];
    Uy=[dUy(i,:) zvec zvec zvec zvec];
    Vx=[zvec dVx(i,:) zvec zvec zvec];
    Vy=[zvec dVy(i,:) zvec zvec zvec];
    Wx=[zvec zvec dWx(i,:) zvec zvec];
    Wy=[zvec zvec dWy(i,:) zvec zvec];
    Px=[zvec zvec zvec dPx(i,:) zvec];
    Py=[zvec zvec zvec dPy(i,:) zvec];
    Rx=[zvec zvec zvec zvec dRx(i,:)];
    Ry=[zvec zvec zvec zvec dRy(i,:)];
    p =[zvec zvec zvec P(i,:)   zvec];
    r =[zvec zvec zvec zvec   R(i,:)];
    
    % Strain transformation matrix [T]: Eq(11)
    C=[Ux; Uy; Vx; Vy; Wx; Wy; Px; Py; Rx; Ry; p; r];

    % Integrated thickness terms
    h=sum(lamDATA.thk);
    if(lamDATA.tDeg > 0)
        h=thk_legendre(lamDATA.tDeg,psiI(i),etaI(i))*lamDATA.thk_coeff'; % assume symmetric distribution
        lamDATA.thk(1:nply)=h/nply;  % assumes even ply %s
    end
    
    % Z offset of element
    z0=zI(i);
    
    % Limit calls to constant psi for speed
    %if( psiI(i) ~= psi_current )
        psi_current = psiI(i);
        % Laminate ABD-Matrix: Eq(11)
        abd=get_abd(panXY, psiI(i), etaI(i), z0, lamDATA, smearKM);

        A=abd(1:3, 1:3);
        B=abd(1:3, 4:6);
        D=abd(4:6, 4:6);
        
        % Plate A-Matrix including transverse shear stiffness
        A=[ A(1,1) A(1,2) A(1,3)  0    0;
            A(2,1) A(2,2) A(2,3)  0    0;
            A(3,1) A(3,2) A(3,3)  0    0;
            0      0      0     Kshear*A(3,3)    0;
            0      0      0       0   Kshear*A(3,3)];
        
        % Plate B-Matrix including transverse shear stiffness
        B=[ B(1,1) B(1,2) B(1,3)  0    0;
            B(2,1) B(2,2) B(2,3)  0    0;
            B(3,1) B(3,2) B(3,3)  0    0;
            0      0      0      B(3,3)    0;
            0      0      0       0      B(3,3)];        
        
        % Plate D-Matrix including transverse shear stiffness
        D=[ D(1,1) D(1,2) D(1,3)  0    0;
            D(2,1) D(2,2) D(2,3)  0    0;
            D(3,1) D(3,2) D(3,3)  0    0;
            0      0      0         D(3,3)   0;
            0      0      0       0        D(3,3)];
        
%             0      0      0     0*D(3,3)   0;
%             0      0      0       0   0*D(3,3)];        
    %end
    
    % Jacobian
    J=jacob2D_Iso( psiI(i), etaI(i), panXY );
    JI=inv(J);
    j11=JI(1,1); j12=JI(1,2); j21=JI(2,1); j22=JI(2,2);

    % Strain transformation matrix
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
    
    % Isotropic material/smeared stiffness
%     A0 = C'*T0'*A*T0*C*h;
%     A1 = ( C'*T0'*A*T1*C + C'*T1'*A*T0*C )*h^2/2;
%     A2 = ( C'*T1'*A*T1*C )*h^3/12;

    % Orthotropic Material
    A0 = C'*T0'*A*T0*C;
    A1 = ( C'*T0'*B*T1*C + C'*T1'*B*T0*C );  % Ignore for symm lam
    % replace A in A1 term with B-matrix later
    %A1 = 0;
    A2 = ( C'*T1'*D*T1*C );
            
    % Stiffness matrix
    K = K + (A0 + A1 + A2)*det(J)*gl_wt*1.02155;

    % Displacement basis functions
    u=[U(i,:) zvec zvec zvec zvec];
    v=[zvec V(i,:) zvec zvec zvec];
    w=[zvec zvec W(i,:) zvec zvec];
    p=[zvec zvec zvec P(i,:) zvec];
    r=[zvec zvec zvec zvec R(i,:)];
       
    H = [u; v; w; p; r];
    
    Z0=[1 0 0 0 0;
        0 1 0 0 0;
        0 0 1 0 0];
    
    Z1=[0 0 0 1 0;
        0 0 0 0 1;
        0 0 0 0 0];
    
    %A0 = H'*Z0'*Z0*H*h; % orig 
    %A0 = H'*Z0'*Z0*H*abs(z+h);
    %A0 = H'*Z0'*Z0*H*(z+h);

    % A1 = H'*(Z0'*Z1+ Z1'*Z0)*H*h^2/4; % ignore term for symm lam

    %A2 = H'*Z1'*Z1*H*h^3/12;
    %A2 = H'*Z1'*Z1*H*abs(z+h)^3/12;
    %A2 = H'*Z1'*Z1*H*(z0^3-(z0-h)^3)/3; %orig 
    
    % Mass Matrix: Eq(24)
    %M = M + rhoh*(A0 + 0*A1 + A2)*det(J)*gl_wt;
    %M = M + rho*(A0 + 0*A1 + A2)*det(J)*gl_wt; %orig 
    
    % Integrate M through section
    z=zeros(nply);
    z(1) = -h/2.;
    for k=1:nply
        z(k+1) = z(k)+lamDATA.thk(k);
    end
    z=z+z0;
    
%     A0=zeros(ndof,ndof);
%     A2=zeros(ndof,ndof);
%     
%     for k=1:nply
%         A0 = A0 + H'*(Z0.'*Z0)*H*(z(k+1) - z(k));
%         A2 = A2 + H'*(Z1.'*Z1)*H*(z(k+1)^3-z(k)^3)/3;
%     end
%     M = M + rho*(A0 + 0*A1 + A2)*det(J)*gl_wt;

    A0=zeros(ndof,ndof);
    A1=zeros(ndof,ndof);
    A2=zeros(ndof,ndof);
    
    for k=1:nply
        rho=lamDATA.rho_ply(k)/(12.*32.2);
        if smearKM==0
            A0 = A0 + rho*(H'*Z0')*(Z0*H)*(z(k+1) - z(k));
            %A1 = A1 + rho*((H'*Z0')*(Z1*H)+(H'*Z1')*(Z0*H))*(z(k+1) - z(k));
            A2 = A2 + rho*(H'*Z1')*(Z1*H)*(z(k+1)^3-z(k)^3)/3;
        else
            A0 = A0 + rho*(H'*Z0')*(Z0*H)*lamDATA.thk(k);
            A2 = A2 + rho*(H'*Z1')*(Z1*H)*((z0+lamDATA.thk(k)/2)^3-(z0-lamDATA.thk(k)/2)^3)/3.;
        end
    end
    M = M + (A0 + A1 + A2)*det(J)*gl_wt*1.02155;

    % Product of density*ply_thk for all plies
    rhoh=(lamDATA.rho_ply*lamDATA.thk.');
    Wt = Wt + rhoh*det(J)*gl_wt*1.02155;
    A=A+det(J)*gl_wt;

    % Integrate M through section - using ZZ matrix 
%    ZZt=zeros(5,5);
    %Mp=zeros(ndof,ndof);
%     for k=1:nply
%         ZZ=[z(k+1)-z(k)       0                0       (z(k+1)^2)/2-(z(k)^2)/2            0;
%                 0           z(k+1)-z(k)        0               0                (z(k+1)^2)/2-(z(k)^2)/2;
%                 0               0        z(k+1)-z(k)           0                          0;
%     (z(k+1)^2)/2-(z(k)^2)/2     0              0       (z(k+1)^3)/3-(z(k)^3)/3            0;
%             0        (z(k+1)^2)/2-(z(k)^2)/2   0               0               (z(k+1)^3)/3-(z(k)^3)/3];
%     %end

        %M = M + lamDATA.rho_ply(k)/(12.*32.2)*H'*ZZ*H*det(J)*gl_wt; 
        %Mp = Mp + lamDATA.rho_ply(k)/(12.*32.2)*(H.'*ZZ*H); 
       
        %ZZt = ZZt+ZZ*lamDATA.rho_ply(k)/(12.*32.2);
        %rhoh = rhoh+lamDATA.rho_ply(k)*lamDATA.thk(k);
        
        
end
%M = M + Mp*det(J)*gl_wt;
%M  = M + H'*ZZt*H*det(J)*gl_wt;