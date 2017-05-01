function [xyzI]=getWingZCoords(Isurf,psiI,etaI,panXY,flip_normal,plot_flag)
%getWingZCoords returns z-coordinates for integration mesh (psiI, etaI)
%   Detailed explanation goes here

addpath('C:\Users\hensomc\Documents\MATLAB\geom3d\geom3d');
xyzI=zeros(length(psiI),3);

zmax=max(max(Isurf.vertices(Isurf.faces(:,:),3)));
zmin=min(min(Isurf.vertices(Isurf.faces(:,:),3)));
ztol=0.010;

%xyzI=zeros(length(Isurf.vertices),3);

% Loop over integration points
for i=1:length(psiI)
%parfor i=1:length(psiI)
    [xI,yI]=ISO_st_to_xy(psiI(i),etaI(i),panXY);

    plane=[];
    c=[];
    
    for j=1:length(Isurf.faces)
        % find face nearest xI,yI
        in = pointInPoly([xI,yI],Isurf.vertices(Isurf.faces(j,:),1:2));
        if(in == 1)
            %disp('point I found in polygon J:'); i,j
            
            plane=createPlane(Isurf.vertices(Isurf.faces(j,1:3),:));
            if(flip_normal >= 1)
                plane=[plane(1:3) -plane(4:6) plane(7:9)];
            end
            
            % alt method: use plane centroid & input normals - needs work
%             c = polygonCentroid3d(Isurf.vertices(Isurf.faces(j,1:4),:));
%             normal = [Isurf.nx(j) Isurf.ny(j) Isurf.nz(j)];
%             plane=createPlane(c, normal);
%             plot3(c(1),c(2),c(3),'o','MarkerSize',12,'MarkerFaceColor','r')
%             quiver3(c(1),c(2),c(3),normal(1),normal(2),normal(3))
            
            % point on plane
            %xyzI(i,1:3) = planePoint(plane, [xI,yI]);
            
            
            % alt method, intersect a line with plane           
            if flip_normal >=1
                z1=-10000;
            else
                z1=10000;
            end
            normal = [Isurf.nx(j) Isurf.ny(j) Isurf.nz(j)];
            c = polygonCentroid3d(Isurf.vertices(Isurf.faces(j,1:4),:));          
            %[I,check] = plane_line_intersect(normal,c,[xI yI 0], [xI yI z1]);
            [I,check] = plane_line_intersect(normal,c,[xI yI -10000], [xI yI 1000]);
            if check==1
                %                 if abs(I(3)) < .0001
                %                     I=[I(1) I(2) 0];
                %                 end
                %xyzI(i,:)=I;
                xyzI(i,1:3)=[xI yI I(1,3)];
            end

            break
        end
    end
    
    % future - investigate a different approach
    %[projection]=[Points]*null([normal vector])
    
    % Pt not found in polygons
    if j==length(Isurf.faces)
        %disp('point not mapped, z-value set to 0.0');xI,yI
        if plot_flag>0
            disp('point not mapped, z-value set to previous projection');xI,yI
        end
        xyzI(i,:)=[xI yI xyzI(i-1,3)];
    end
    
    % Reset bad projections to x,y,0.0
    if xyzI(i,3)>zmax+ztol || xyzI(i,3)<zmin-ztol
        if plot_flag>0
            disp('Bad projection at : , z-value set to 0.0'); xyzI(i,:)
        end
        xyzI(i,:)=[xI yI 0.0];
    end
    
end
%xyzI

if plot_flag>0
    %scatter3(xyzI(:,1)',xyzI(:,2)',xyzI(:,3)','o')
    scatter3(xyzI(:,1)',xyzI(:,2)',xyzI(:,3)','o','MarkerFaceColor','b');
    %figure,plot3(x,y,z,'o')
end