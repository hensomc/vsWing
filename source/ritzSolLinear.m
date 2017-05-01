function [C, ritzDisp, ritzEps]=ritzSolLinear( solDATA, panXY, lamDATA, psiBCMesh, etaBCMesh, surfI, Kq, Fq, T, plot_flag )

% Message
if(plot_flag >=1)
    disp('Solving for displacements...');
end

% Recover parameters
pDeg=solDATA.pDeg;
ipoltype=solDATA.ipoltype;


% Solve
q=Kq\Fq;
C=T*q;

% Ritz displacement vector - Use BCMesh to compare Ritz sol with FEM mesh
[F,Fx,Fy,Fxy,Fxx,Fyy] = FdF(ipoltype, pDeg, psiBCMesh(:),etaBCMesh(:));
dim=(pDeg+1)^2;

Uout=F*C(1:dim);
Vout=F*C(dim+1:2*dim);
Wout=F*C(2*dim+1:3*dim);
Pout=F*C(3*dim+1:4*dim);
Rout=F*C(4*dim+1:5*dim);

ritzDisp=[Uout Vout Wout Pout Rout];  % Dim[ (pDeg+1)^2, 5]

% compute strains - will need Jacobian:
ritzEps=getRitzStrain( solDATA, panXY, lamDATA, surfI, C, F,Fx,Fy,Fxy,Fxx,Fyy, psiBCMesh(:),etaBCMesh(:) );
