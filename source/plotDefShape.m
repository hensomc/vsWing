function plotDefShape(panXY,solDATA,lamDATA,psiIMesh,etaIMesh,C,plot_flag)

% Recover data from structures
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);


% Get titles
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);

%% plot deformed shape
ncases = 1;

% One figure for all deflected shapes
if( plot_flag == 0 )
    return
else
    %figure
end

% Get mesh in x,y space
[xIMesh,yIMesh] = meshTransform(panXY, psiIMesh, etaIMesh);

% Deformed shape for each case
for i=1:ncases
    
    % Compute displacement components
    [FOUT] = FdF(ipoltype, pDeg, psiIMesh(:), etaIMesh(:));
    dim=(pDeg+1)^2;
    uout=FOUT*C(1:dim);
    vout=FOUT*C(dim+1:2*dim);
    wout=FOUT*C(2*dim+1:3*dim);
    %disp('size of C=');size(C)
    %disp('size of Wout=');size(Wout)
    

    %% Displaced shape using points
    figure
    %plot(xDMesh(:),yDMesh(:),'r.','linewidth',3);
    plot(xIMesh(:)+uout,yIMesh(:)+vout,'r.','linewidth',3);
    %daspect([1 1 1]);
    %axis equal; axis tight;
    title({'\bfInplane Deformations'});
    hold off
    
    %% Deformed shape with filled contour plot
    %resDisp=sqrt(uout'*uout +vout'*vout)
    uu=reshape(uout,size(xIMesh));
    vv=reshape(vout,size(xIMesh));
    ww=reshape(wout,size(xIMesh));
    
    % Inplane displacements
    figure
    %res=sqrt(power(uu,2)+power(vv,2));
    surfc(xIMesh,yIMesh,vv);  % inplane loads
    %daspect([1 1 1]);
    view(2)
    title({'\bfRitz Analysis Deformed Shape';'v-Displacement Contours'},'fontsize',14);
    hold off
    
    % 10X Scaled inplane displacemnts
    figure
    surf(xIMesh+25*uu,yIMesh+25*vv,vv);
    title({'\bfRitz Analysis 10X Deformed Shape';'v-Displacement Contours'},'fontsize',14);
    view(2)
    hold off

    % Transverse displacements   
    figure
    surfc(xIMesh,yIMesh,ww);  % transverse loads   
    %daspect([abs(min(panXY(:,1))-max(panXY(:,1))) abs(min(panXY(:,2))-max(panXY(:,2))) 1]);
    title({'\bfRitz Analysis Deformed Shape';'w-Displacement Contours'},'fontsize',14);
    hold off
    
    return

    % Plot shapes
    if( plot_flag == 1 || plot_flag == 3 )
        %subplot(2,3,i)
        WW=reshape(Wout,size(xIMesh));
        %disp('size of WW='); size(WW)
        %figure
        surfc(xIMesh,yIMesh,WW)
        shading interp
        colorbar
        title1='\bfRitz Analysis: w Deflected Shape';
        caseTitle=strcat('Case=',int2str(i));
        title( {title1;panTitle;caseTitle} );
        %title({['\bfCase ',int2str(i)]; ['Case=',num2str(i,'%d')]})
    end
    
    % Plot animated shapes
    if( plot_flag == 2 )
        WW=reshape(Wout,size(xIMesh));
        %disp('size of WW='); size(WW)
        figure
        h=surf(xIMesh,yIMesh,WW);
        shading interp
        %axis([min(min(xIMesh)) max(max(xIMesh)) min(min(yIMesh)) max(max(yIMesh)) min(min(WW)) max(max(WW))]);
        axis([min(min(xIMesh)) max(max(xIMesh)) min(min(yIMesh)) max(max(yIMesh)) min(min(WW)) max(max(WW))]); % in-plane displ
        %colorbar
        
        % animate w-displacement
        az=-37.5;
        for iscale=0 : 0.05 : 1.0;
            set(h, 'xdata', xIMesh, 'ydata', yIMesh, 'zdata', iscale*WW)
            view(az, 30)
            drawnow
            pause(0.05)
            az=az+3;
        end
        for iscale=1 : -0.05 : -1.0;
            set(h, 'xdata', xIMesh, 'ydata', yIMesh, 'zdata', iscale*WW)
            view(az, 30)
            drawnow
            pause(0.05)
            az=az+3;
        end
        for iscale=-1 : 0.05 : 1.0;
            set(h, 'xdata', xIMesh, 'ydata', yIMesh, 'zdata', iscale*WW)
            view(az, 30)
            drawnow
            pause(0.05)
            az=az+3;
        end
        % Modify azimuth (horizontal rotation) and update drawing
        %       for az = -37.5 : .5 : 30
        %           view(az, 30)
        %           drawnow
        %       end
    end
end
