function plotFibCurvatureSpace(panXY,Rallow,phiSkew)
% Plots fiber curvature in +x space of rectangular panel of width = a

a=abs(panXY(2,1)-panXY(1,1));
b=abs(panXY(2,2)-panXY(1,2));
x=linspace(0,a/2);
y=linspace(0,b/2);

T0=linspace(0,90);
T1=linspace(-90,90);

% Matrix of curvatures over space T0,T1
for i=1:length(T0)
    for j=1:length(T1)
        % Square plate curvatures
        fibK = getFibCurvature(a,T0(i),T1(j),x);
        fibK_min(i,j) = min(min(abs(fibK)));
        fibK_max(i,j) = max(max(abs(fibK)));

        % Skew plate curvatures
        fibKskew = getFibCurvatureSkew(a,T0(i),T1(j),x,phiSkew);
        fibKskew_min(i,j) = min(min(abs(fibKskew)));
        fibKskew_max(i,j) = max(max(abs(fibKskew)));         
    end
end

% Create mesh for region plot and contour plot
%[xx,yy]=meshgrid(T0,T1);
[xx,yy]=meshgrid(T1,T0);
X=xx(:);Y=yy(:);F=fibK_min(:); 

% For ref only
max(F)  % find max curvature
min(F)  % find min curvature

% Find curvatures less than input allowable
I=find(F<=1/Rallow);

% Get same results for radius of curvature: reset large/singular values to
% Rallow for ease of plotting
fibR_min=1./fibK_min;
R=find(fibR_min >= Rallow);
%fibR_min(R)=Rallow; - use this to improve shading

% Plot allowable fiber curvature design space as region of points
figure; plot(X(I),Y(I),'r.')
title( {['\bfDesign Space for Shifted Curvlinear Fibers']; ['R_{allow}= ',num2str(Rallow)]},'fontsize',14 )
xlabel('\bfT_1','fontsize',12)
ylabel('\bfT_0','fontsize',12)

% Plot design space as contour
figure; contourf(xx,yy,fibK_min); colorbar
title( {['\bfContour of Shifted Fiber Curvatures']; ['R_{allow}= ',num2str(Rallow)]},'fontsize',14 )
xlabel('\bfT_1','fontsize',12)
ylabel('\bfT_0','fontsize',12)

% Plot inverse of curvature = radius of curvature
h1=figure; cc=contourf(xx,yy,reshape(fibR_min,size(xx))); colorbar
%clabel(cc);
title( {['\bfDesign Space for Shifted Curvlinear Fibers']; ['R_{allow}= ',num2str(Rallow)]},'fontsize',14 )
xlabel('\bfT_1','fontsize',12)
ylabel('\bfT_0','fontsize',12)

% Radius of curvature - grayscale
h2=figure; cc=contourf(xx,yy,reshape(fibR_min,size(xx))); colorbar
title( {['\bfDesign Space for Shifted Curvlinear Fibers']; ['R_{allow}= ',num2str(Rallow)]},'fontsize',14 )
xlabel('\bfT_1','fontsize',12)
ylabel('\bfT_0','fontsize',12)
set(gcf,'colormap',gray);

% =============== SKEW RESULTS ==================
F2=fibKskew_min(:);

% Find curvatures less than input allowable
I2=find(F2<=1/Rallow);
% Plot allowable fiber curvature design space as region of points
figure; plot(X(I2),Y(I2),'r.')
title( {['\bfDesign Space for Shifted Curvlinear Fibers']; ['R_{allow}= ',num2str(Rallow)]},'fontsize',14 )
xlabel('\bfT_1','fontsize',12)
ylabel('\bfT_0','fontsize',12)
