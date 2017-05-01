function [surfI]=getDegenSurf_fvc(degenGeom, isurf, plot_flag)

% Global history figure handle - need to preserve
global h_histFig;
global PlotThkSurf;

%wingSurface recovers 3D face data for upper and lower wing surfaces
%   Detailed explanation goes here

% usurf.zI=zeros(1,length(psiI));
% lsurf.zi=zeros(1,length(psiI));

if(size(degenGeom) <= 0)
    return
end

% surface coordinates
x=12*degenGeom(1).surf.x;
y=12*degenGeom(1).surf.y;
z=12*degenGeom(1).surf.z;

% normals
nx=degenGeom(1).surf.nx;
ny=degenGeom(1).surf.ny;
nz=degenGeom(1).surf.nz;

% Index of lower/upper surface intersection
s=size(x);
s1=round(s(2)/2);
s2=s(2);

% Upper surface 
xu=x(:,s1:s2);
yu=y(:,s1:s2);
zu=z(:,s1:s2);

% Draw upper surface
if isurf == 1
    %if plot_flag>0
%         figure;
        hu=surf(xu,yu,zu);
        %hold on;
%         hold off;
        usurf=surf2patch(hu);
        %if plot_flag == 0
            if gcf==h_histFig
                set([hu],'Visible','off');
            elseif ismember(gcf,PlotThkSurf)==1
                set([hu],'Visible','off');
            else
                set([hu,gcf],'Visible','off');
            end
        %end
        
%         figure;
%         M=mesh(xu,yu,zu);
%         set(M,'facecolor','none');
%         hold on;
%     else
%         %usurf=surf2patch(surf(xu,yu,zu));
%         s=surf(xu,yu,zu);
%         usurf=surf2patch(s);
%         %set(usurf,'Visible','off');
%         if gcf==h_histFig
%         else
%             set([s,gcf],'Visible','off');
%         end
% 
%     end
end
% normals
n=size(nx);
n1=round(n(2)/2)+1;
n2=n(2);

usurf.nx=reshape(nx(:,n1:n2),1,numel(nx(:,n1:n2)));
usurf.ny=reshape(ny(:,n1:n2),1,numel(ny(:,n1:n2)));
usurf.nz=reshape(nz(:,n1:n2),1,numel(nz(:,n1:n2)));

% Upper surface integration points
% [xyzIu]=getWingZCoords(usurf,psiI,etaI,panXY,0);
% plot3(xyzIu(:,1)',xyzIu(:,2)',xyzIu(:,3)','ro','markersize',12)
% usurf.zI=xyzIu(:,3)';
% hold off;

% Lower surface
xl=x(:,1:s1);
yl=y(:,1:s1);
zl=z(:,1:s1);

if isurf ~= 1
%     if plot_flag>0
%        figure;
        hl=surf(xl,yl,zl);
%        hold on;
%        hold off;
        lsurf=surf2patch(hl);
        %if plot_flag == 0
        if gcf==h_histFig
            set([hl],'Visible','off');
        elseif ismember(gcf,PlotThkSurf)==1
            set([hl],'Visible','off');
        else
            set([hl,gcf],'Visible','off');
        end
        %end
        
%        figure
%         M=mesh(xl,yl,zl);
%         set(M,'facecolor','none');
%         hold on;
%     else
%         %lsurf=surf2patch(surf(xl,yl,zl));
%         s=surf(xu,yu,zu);
%         lsurf=surf2patch(s);
%         if gcf==h_histFig;
%         else
%             set([s,gcf],'Visible','off');
%         end
%     end

end

lsurf.nx=reshape(nx(:,1:n1),1,numel(nx(:,1:n1)));
lsurf.ny=reshape(ny(:,1:n1),1,numel(nx(:,1:n1)));
lsurf.nz=reshape(nz(:,1:n1),1,numel(nx(:,1:n1)));

% Lower surface integration points
% [xyzIl]=getWingZCoords(lsurf,psiI,etaI,panXY,1);
% plot3(xyzIl(:,1)',xyzIl(:,2)',xyzIl(:,3)','ro','markersize',12)
% lsurf.zI=xyzIl(:,3)';
% hold off;

if isurf==1
    surfI=usurf;
else
    surfI=lsurf;
end
