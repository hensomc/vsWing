function [ layer_thk ] = getLayerThk(tDeg,psiI,etaI,lamDATA,plot_flag,titleIn)
%getLayerThk - returns layer_thk for surface with polynomial thickness
%              distribution - uses Legendre polynomial

% Default return argument
layer_thk=[];

% Return if not variable thickness layer
if tDeg <= 0
    return
end

% dim1=length(psiI) ~ (Mg+1)^2
% dim2=(tDeg+1)^2
% layer_thk = [A]*{thk_coeff}; A(dim1,dim2)*thk_coeff{dim2,1) 
layer_thk = thk_legendre(tDeg,psiI,etaI)*lamDATA.thk_coeff';

plotLayerThk(psiI,etaI,layer_thk, titleIn)

function plotLayerThk(psiI,etaI,layer_thk, titleIn)
% Plot the thickness as a surface function

% Global variable to capture multiple thickness surface plots
global PlotThkSurf;
global iPlotThkSurf;

%if( plot_flag == 1 )
    dim=sqrt(length(psiI));
    lay_thk_reshape=reshape(layer_thk,dim,dim);
    psiI2=reshape(psiI,dim,dim);
    etaI2=reshape(etaI,dim,dim);
    %disp('size lay_thk = '); size(lay_thk_reshape)
    h=figure;
    PlotThkSurf(iPlotThkSurf)=h;
    h=surfc(psiI2,etaI2,lay_thk_reshape);
    shading interp
    %colorbar
    
    % ----------------------------
    % add transparent plane at z=0
    % ----------------------------
    z_threshold = 0.0;
    
    % Obtain the limits of the axes
    yp = get(gca,'Ylim');
    xp = get(gca,'Xlim');
    
    % Use the axes x and Y limits to find the co-ordinates for the patch
    x1 = [ xp(1) xp(2) xp(2) xp(1)];
    y1 = [ yp(1) yp(1) yp(2) yp(2)];
    
    % Patch is parallel to XY plane, Z coordinate = z_threshold
    z1 = ones(1,numel(x1))* z_threshold;  % creates a 1x4 vector representing the Z coordinate values 
    p = patch(x1,y1,z1, 'b');
    
    % Set the Face and edge transparency to 0.2 using the following properties
    set(p,'facealpha',0.2)
    set(p,'edgealpha',0.2)
    
    if length(titleIn)>0
        title(titleIn);
    else
        title('\bf Thickness Distribution')
    end
    hold off
    %pause
    
%end

