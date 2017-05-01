% Stiffness and Mass Matrix for Web with Transverse Shear Terms
function [K,M,Wt]=KMweb(pDeg,intMethod,panXY,lamDATA,psiI,etaI,zIU,zIL,gl_wts,B,eqW,plot_flag)

if(plot_flag>0)
    disp('Calculating stiffness matrix...');
end
ndof=5*(pDeg+1)^2;
K=zeros(ndof,ndof);M=zeros(ndof,ndof);
zvec=zeros(1,(pDeg+1)^2);

Wt=0.0;

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

% Trot = property transform based  on beam x-dir
ang = 90*pi/180.
c=cos(ang);
s=sin(ang);

c2=c*c;
s2=s*s;
c3=cos(2*ang);
s3=sin(2*ang);
s4=0.5*s3;

Trot=[ c2  s2  s3  0   0;
       s2  c2 -s3  0   0;
      -s4  s4  c3  0   0;
        0   0   0  c  -s;
        0   0   0  s   c];

% inv(Trot) = Trotp = Trot(-theta);
c=cos(-ang);
s=sin(-ang);

c2=c*c;
s2=s*s;
c3=cos(2*ang);
s3=sin(2*ang);
s4=0.5*s3;

Trotp=[ c2  s2  s3  0   0;
        s2  c2 -s3  0   0;
       -s4  s4  c3  0   0;
         0   0   0  c  -s;
         0   0   0  s   c];

nz = 10;  % assume 10 subdivisions of web

for i=1:length(psiI);
    dz=(zIU(i)-zIL(i))/nz;
    z=zIL(i)+dz/2;
    zz=zIL(i);

    for j=1:nz-1    
        gl_wt=1.0;
        if(intMethod == 2)  %gaussian quadrature weights (in 2D)
            gl_wt=gl_wts(i);
        end
        
        % Strain basis functions
        Ux=[dWx(i,:) zvec zvec zvec zvec];
        Uy=[dWy(i,:) zvec zvec zvec zvec];
        Vx=[zvec dVx(i,:) zvec zvec zvec];
        Vy=[zvec dVy(i,:) zvec zvec zvec];
        Wx=[zvec zvec dUx(i,:) zvec zvec];
        Wy=[zvec zvec dUy(i,:) zvec zvec];
        Px=[zvec zvec zvec zvec zvec];
        Py=[zvec zvec zvec zvec zvec];
        Rx=[zvec zvec zvec zvec zvec];
        Ry=[zvec zvec zvec zvec zvec];
        p =[zvec zvec zvec zvec   zvec];
        r =[zvec zvec zvec zvec   zvec];
        
        % Strain transformation matrix [T]: Eq(11)
        C=[Ux; Uy; Vx; Vy; Wx; Wy; Px; Py; Rx; Ry; p; r];
        
        % Integrated thickness terms
        h=sum(lamDATA.thk);
        
        % Limit calls to constant psi for speed
        %if( psiI(i) ~= psi_current ) - requires properties const with psi
        psi_current = psiI(i);
        % Laminate ABD-Matrix: Eq(11)
        %abd=get_abd(panXY, psiI(i), etaI(i), 0, lamDATA);
        abd=get_abd(panXY, psiI(i), etaI(i), z, lamDATA);
        
        A=abd(1:3, 1:3);
        D=abd(4:6, 4:6);
        
        %Plate A-Matrix including transverse shear stiffness
%                 A=[ A(1,1) A(1,2) A(1,3)  0    0;
%                     A(2,1) A(2,2) A(2,3)  0    0;
%                     A(3,1) A(3,2) A(3,3)  0    0;
%                     0      0      0      A(3,3)    0;
%                     0      0      0       0   A(3,3)];
%         A=[ A(1,1) 0      0       0      0;
%             0      0      0       0      0;
%             0      0   A(3,3)     0      0;
%             0      0      0     A(3,3)   0;
%             0      0      0       0   A(3,3)];
        A=[ A(1,1) 0      0       0      0;
            0      0      0       0      0;
            0      0      0       0      0;
            0      0      0       0      0;
            0      0      0       0   A(3,3)];
        
        A=Trot*A*Trotp;
        
        % Plate D-Matrix including transverse shear stiffness
        %         D=[ D(1,1) D(1,2) D(1,3)  0    0;
        %             D(2,1) D(2,2) D(2,3)  0    0;
        %             D(3,1) D(3,2) D(3,3)  0    0;
        %             0      0      0       D(3,3)   0;
        %             0      0      0       0      D(3,3)];
%         D=[ D(1,1) 0      0       0      0;
%             0      0      0       0      0;
%             0      0    D(3,3)     0      0;
%             0      0      0   D(3,3)     0;
%             0      0      0       0      D(3,3)];
        D=[ 0      0      0       0      0;
            0      0      0       0      0;
            0      0      0       0      0;
            0      0      0       0     0;
            0      0      0       0      0];        D=Trot*D*Trotp;
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
        % A1 = ( C'*T0'*A*T1*C + C'*T1'*A*T0*C )*h^2/2;  % Ignore for symm lam
        % replace A in A1 term with B-matrix later
        A1 = 0;
        A2 = ( C'*T1'*D*T1*C );
        
        % Beam equivalent width (normalized)
        if eqW(3) == 0 % spar: beam axis || to psi
            bmW=0.5*eqW(1)*(1-etaI(i))+0.5*eqW(2)*(1+etaI(i));
        elseif eqW(3) == 1 % rib: beam axis || to eta
            bmW=0.5*eqW(1)*(1-psiI(i))+0.5*eqW(2)*(1+psiI(i));
        end
        
        % Stiffness matrix
        K = K + bmW*(A0 + 0*A1 + 0*A2)*det(J)*gl_wt;
        
        % Density
        rhoh=(lamDATA.rho_ply*lamDATA.thk.')/(12.*32.2);
        rho=lamDATA.rho_ply(1)/(12.*32.2);
        
        % Displacement basis functions
        u=[U(i,:) zvec zvec zvec zvec];
        v=[zvec V(i,:) zvec zvec zvec];
        w=[zvec zvec W(i,:) zvec zvec];
        p=[zvec zvec zvec P(i,:) zvec];
        r=[zvec zvec zvec zvec R(i,:)];
        
        H = [w; v; u; p; r];
        %H = [v; w; u; p; r];
        
%         Z0=[1 0 0 0 0;
%             0 1 0 0 0;
%             0 0 1 0 0];
%         
%         Z1=[0 0 0 1 0;
%             0 0 0 0 1;
%             0 0 0 0 0];
%         
%         %A0 = H'*Z0'*Z0*H*h;
%         A0 = H'*Z0'*Z0*H*abs(z);
%         % A1 = H'*(Z0'*Z1+ Z1'*Z0)*H*h^2/4; % ignore term for symm lam
%         A1 = 0;
%         %A2 = H'*Z1'*Z1*H*h^3/12;
%         A2 = H'*Z1'*Z1*H*abs(z)^3/12;
%         
%         % Mass Matrix: Eq(24)
%         %M = M + rhoh*(A0 + 0*A1 + A2)*det(J)*gl_wt;
%         M = M + rho*(A0 + 0*A1 + A2)*det(J)*gl_wt;


        ZZ=[z-zz       0                0       (z^2)/2-(zz^2)/2            0;
                0           z-zz        0               0                (z^2)/2-(zz^2)/2;
                0               0        z-zz           0                          0;
    (z^2)/2-(zz^2)/2     0              0       (z^3)/3-(zz^3)/3            0;
            0        (z^2)/2-(zz^2)/2   0               0               (z^3)/3-(zz^3)/3];

        M = M + rho*bmW*H'*ZZ*H*det(J)*gl_wt;  
        
        Wt = Wt + rho*(z-zz)*det(J)*gl_wt;
        zz=z;
        z=z+dz;
        
    end
end
