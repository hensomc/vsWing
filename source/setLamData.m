function [ lamDATA ] = setLamData( ilam, tDeg, thk_coeff )
%function [ lamDATA ] = setLamData( ilam, pThk )
%setLamData Sets the laminate data structure
%   ilam = input laminate ID
%   tDeg = degree of legendre thickness polynomial
%   Polynomial thickness coefficients [c0, c1, c2, c3, c4]; t=c0*t0 +c1*x +c2*y + c3*xy


%   lamDATA = returned data structure
%       lamDATA.imat      = ply material IDs;
%       lamDATA.theta0    = ply theta0;
%       lamDATA.theta1    = ply theta1;
%       lamDATA.thk       = ply thickneses ;
%       lamDATA.mat_names = material names;  
%       lamDATA.phi_rot   = ply rotation phi;

%   material data
%       lamDATA.e1   = longitudinal modulus;
%       lamDATA.e2   = transverse modulus;
%       lamDATA.nu12 = poissons ratio;
%       lamDATA.g12  = shear modulus ;
%       lamDATA.rho   = density;  

% Laminate definition
[imat, theta0, theta1, thk, mat_names, phi_rot] = readLam(ilam);

lamDATA.ilam=ilam;
lamDATA.imat=imat;
lamDATA.phi_rot=phi_rot;

lamDATA.theta0=theta0;
lamDATA.theta1=theta1;
lamDATA.thk=thk;
lamDATA.mat_names=mat_names;    

% Laminae material properties
[e1,e2,nu12,g12,rho,mtype] = get_mat_props(mat_names);
lamDATA.e1   = e1;
lamDATA.e2   = e2;
lamDATA.nu12 = nu12;
lamDATA.g12  = g12 ;
lamDATA.rho  = rho;
lamDATA.mtype= mtype;

for iply=1:length(thk)
    rho_ply(iply)=rho(imat(iply));
end;
lamDATA.rho_ply=rho_ply;

lamDATA.tDeg = tDeg;

lamDATA.thk_coeff = thk_coeff;

% Thickness coeff's for upper/lower surface layers
% lamDATA.usurf_thk_coeff = zeros(4,4);
% lamDATA.lsurf_thk_coeff = zeros(4,4);
% lamDATA.usurf_thk_coeff(1,:) = thk_coeff;
% lamDATA.lsurf_thk_coeff(1,:) = thk_coeff;
lamDATA.usurf_thk_coeff = thk_coeff;
lamDATA.lsurf_thk_coeff = thk_coeff;


% Thk coefficients
% lamDATA.pThk = pThk;
% switch pThk
%     case 0
%         coeff=[1 0 0 0]; % t=c1*t0
%     case 1
%         coeff=[1 1 1 1]; % t=c1*t0 + c2*x + c3*y +c4*x*y
%     case 2
%         coeff=[1 1 1 1 1 1]; % t=c1*t0 + c2*x + c3*y +c4*x*y +c5*x^2 +c5*y^2
% end
% for iply=1:length(pThk)
%     thk_coeff(iply,:)=coeff;
% end

