function preSolPlots(plot_flag, wingDATA, lamDATA, solDATA, psiBCMesh, etaBCMesh, xBCMesh, yBCMesh, psiI, etaI)

if( plot_flag <= 0 )
    return
end

% Recover panel data
panXY = wingDATA.panXY;

% Draw wing geom
drawWingGeom( plot_flag, wingDATA );

% Draw fiber paths and ply courses
refPath=1;  % Draw a reference fiber path for each ply
panFibAxis=1;  % Panel fiber axis variation: ( 0=f(x)|1=f(y) )
draw_panel_fib(plot_flag, panXY, lamDATA, 10, refPath, panFibAxis);  % NumFibPts=10 for speed

% Draw fiber path field
draw_ply_theta_field(panXY,lamDATA,10,refPath,panFibAxis,xBCMesh,yBCMesh);


% Draw fiber curvature space
%plotFibCurvatureSpace(panXY,Rallow,phiSkew);

% Calc & plot lam engr constants as f(x)
%plot_lam_props(solDATA,panXY, lamDATA);

% Contour plot of lam egr constants
%cplot_lam_props(solDATA,panXY, psiBCMesh, etaBCMesh, lamDATA);
%pause

% Draw layer thickness profile as a surface
layer_thk = getLayerThk(lamDATA.tDeg,psiI,etaI,lamDATA,plot_flag,'');

% Plot BC grids
drawBC_grids(plot_flag,psiBCMesh,etaBCMesh,xBCMesh,yBCMesh);

% Plot integration element corners and integ points at centers
if(plot_flag == 3 && intMethod == 1)
    figure; plot(psiIMesh, etaIMesh, 'ko', 'MarkerSize', 2, 'MarkerFaceColor', [0 0 0]);
    axis equal; axis tight; hold on
    plot(psiI, etaI, 'bs', 'MarkerSize', 3);
    title({'\bfNon-Dimensional Integration points'})
    hold off;
end

if(plot_flag >=1)
    Mg=solDATA.Mg;
    order_1d=[Mg+1,Mg+1];
    order_nd = prod ( order_1d(1:2) );
    gl_pts(1,:)=psiI;
    gl_pts(2,:)=etaI;
    gl_grid_display2( order_1d, gl_pts );
end
