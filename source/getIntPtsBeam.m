function [psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPtsBeam(Mg,intMethod,beamXY,panXY,plot_flag)

% Define a fine mesh & integration subdivisions
nISubDiv = 20;                                    % number of integration subdivisions

psiIDomain = linspace(-1,1,nISubDiv);             % integration domain x-dir
etaIDomain = linspace(-1,1,nISubDiv);             % integration domain y-dir
dpsi = psiIDomain(2)-psiIDomain(1);                     % dx = distance between 2 integration points
deta = etaIDomain(2)-etaIDomain(1);                     % dy = distance between 2 integration points
[psiIMesh,etaIMesh] = meshgrid(psiIDomain,etaIDomain);

% Nondimensional integration points - default locations, redfined
%                                     if gauss quadrature integration
psiI  = psiIMesh(:) + dpsi/2;        % x-midpoint of integration element
etaI  = etaIMesh(:) + deta/2;        % y-midpoint of integration element

IIDpsi  = find(psiI<1);           % filter out points on boundary psi = 1
psiI   = psiI(IIDpsi);
IIDeta  = find(etaI<1);           % filter out points on boundary eta = 1
etaI   = etaI(IIDeta);

IIDpsieta= intersect(IIDpsi,IIDeta);
psiI   = psiI(IIDpsieta);
etaI   = etaI(IIDpsieta);

gl_wts = 0;


%% Gauss Legendre points and weights
% get s,t for beam endpts
[s1,t1]=ISO_XY_to_st_Num(beamXY(1,1),beamXY(1,2),panXY);
[s2,t2]=ISO_XY_to_st_Num(beamXY(2,1),beamXY(2,2),panXY);

% beam vector
v=[s2-s1 t2-t1 0];

% plate vector 1-2
[s1p,t1p]=ISO_XY_to_st_Num(panXY(1,1),panXY(1,2),panXY);
[s2p,t2p]=ISO_XY_to_st_Num(panXY(2,1),panXY(2,2),panXY);
u=[s2p-s1p t2p-t1p 0];

CosTheta = dot(u,v)/(norm(u)*norm(v));
Theta = acosd(CosTheta);

if(intMethod == 2)
    % Uniform mesh
%     psiIDomain = linspace(-1,1,nISubDiv);             % integration domain x-dir
%     etaIDomain = linspace(-1,1,nISubDiv);             % integration domain y-dir
    psiIDomain = linspace(s1,s2,nISubDiv);             % integration domain x-dir
    etaIDomain = linspace(t1,t2,nISubDiv);             % integration domain y-dir
    [psiIMesh,etaIMesh] = meshgrid(psiIDomain,etaIDomain);
    
    % 2D Gauss integration mesh
    order_1d=[Mg+1,Mg+1];
    order_nd = prod ( order_1d(1:2) );
    [gl_pts,gl_wts]=gl_grid( 2, order_1d, order_nd );
    psiI=gl_pts(1,:);
    etaI=gl_pts(2,:);
    
    % Use 1D Gauss pts and weights along beam axis
    [gl_pts_1d,gl_wts_1d]=lgwt(Mg+1,-1,1);
    
    % Assign psi & eta based on beam orientation
    tol = 0.001;
    if abs(Theta-0) < tol
        psiI=gl_pts_1d;
        etaI=s1*ones(size(psiI));
    elseif abs(Theta-90.0) < tol
        etaI=gl_pts_1d;
        psiI=s1*ones(size(etaI));
    end
    gl_wts=gl_wts_1d;
end
