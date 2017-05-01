function [projData]=projectPointsToSurf(Isurf,xydata,flip_normal,plot_flag)
% returns xyz-coordinates for 2D pt projected onto 3D surf
%   Detailed explanation goes here

addpath('C:\Users\hensomc\Documents\MATLAB\geom3d\geom3d');
projData.xyz=zeros(length(xydata),3);

% get bounds on z from surface faces
zmax=max(max(Isurf.vertices(Isurf.faces(:,:),3)));
zmin=min(min(Isurf.vertices(Isurf.faces(:,:),3)));
ztol=0.010;


% Loop over  points
for i=1:length(xydata)
%parfor i=1:length(xydata)

    plane=[];
    c=[];
%     x=xydata(i,2);
%     y=xydata(i,3);
    x=xydata(i,1);
    y=xydata(i,2);
   
    for j=1:length(Isurf.faces)
        % find face nearest x,y
        %on = pointInPoly([x,y],Isurf.vertices(Isurf.faces(j,:),1:2));

        edges=[1 2; 2 3; 3 4; 4 1];

        %[in,on] = pointInPoly([x,y],Isurf.vertices(Isurf.faces(j,:),1:2),edges);
        
        [in,on] = inpolygon(x,y,Isurf.vertices(Isurf.faces(j,:),1),Isurf.vertices(Isurf.faces(j,:),2));  %Matlab Fn
        
        if in == 1 || on == 1
            %disp('point I found in polygon J:'); i,j
%             disp('in,on');in,on
            
            plane=createPlane(Isurf.vertices(Isurf.faces(j,1:3),:));
            if(flip_normal >= 1)
                plane=[plane(1:3) -plane(4:6) plane(7:9)];
            end
%             disp('x,y');x,y
%             disp('vertices'),Isurf.vertices(Isurf.faces(j,:),1:2)
            
            % alt method to define plane: plane centroid & plane normal -
            c = polygonCentroid3d(Isurf.vertices(Isurf.faces(j,1:4),:));
            normal = [Isurf.nx(j) Isurf.ny(j) Isurf.nz(j)];
            % gives wrong projection
%             plane=createPlane(c, normal);
%             plot3(c(1),c(2),c(3),'o','MarkerSize',12,'MarkerFaceColor','r')
%             quiver3(c(1),c(2),c(3),normal(1),normal(2),normal(3))
%             
            
            % planePoint method
            %projData.xyz(i,:) = planePoint(plane, [x,y]);
            
            % pointOnPlane method (gives same as above?)
%             projData.xyz(i,:) = projPointOnPlane([x,y,0], plane);
%             projData.xyz(i,1:3);
            
            % alt method, intersect a line with plane
            if flip_normal >=1
                z1=-10000;
            else
                z1=10000;
            end
            %[I,check] = plane_line_intersect(normal,c,[x y 0], [x y z1]);
            [I,check] = plane_line_intersect(normal,c,[x y -10000], [x y 10000]);
            if check==1
%                 if abs(I(3)) < .0001
%                     I=[I(1) I(2) 0];
%                 end
                %projData.xyz(i,:)=I;
                projData.xyz(i,:)=[x y I(1,3)];
            end
            
            break
        end
    end
    if j==length(Isurf.faces)
        if plot_flag>0
            disp('point not mapped, z-value set to previous point offset');x,y, projData.xyz(i-1,3)
        end
        %projData.xyz(i,:)=[x y 0.0];
        projData.xyz(i,:)=[x y projData.xyz(i-1,3)];
    end
    
    % Reset bad projections to x,y,0.0
    if projData.xyz(i,3)>zmax+ztol || projData.xyz(i,3)<zmin-ztol
        if plot_flag>0
            disp('Bad projection at : , z-value set to 0.0'); projData.xyz(i,:)
        end
        projData.xyz(i,:)=[x y 0.0];
    end
    
end
%projData.xyzI

if plot_flag>0
    scatter3(projData.xyz(:,1)',projData.xyz(:,2)',projData.xyz(:,3)','o')
end
%figure,plot3(x,y,z,'o')