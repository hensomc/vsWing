function plotModes(panXY,ipoltype,pDeg,bcMethod,psiIMesh,etaIMesh,C,EgN,plot_flag)
%% plot mode shapes
NMode=min([6 length(EgN)]);

% One figure for all mode shapes
if( plot_flag == 0 )
    return
else
    h1=figure;
    h2=figure;
end

% Get mesh in x,y space
[xIMesh,yIMesh] = meshTransform(panXY, psiIMesh, etaIMesh);

for i=1:NMode
    
    % Compute mode shape at integration points
    [FOUT] = FdF(ipoltype, pDeg, psiIMesh(:), etaIMesh(:));
    Cm=C{i};  % coeff's for ith mode
    iCoeff=(pDeg+1)^2; % index for each DOF
    uout=FOUT*Cm(1:iCoeff);
    vout=FOUT*Cm(iCoeff+1:2*iCoeff);
    wout=FOUT*Cm(2*iCoeff+1:3*iCoeff);
    
    % Plot 2-D modes shapes
    if( plot_flag == 1 || plot_flag == 3 )
        set(0,'CurrentFigure',h1);
        subplot(2,3,i)
        WW=reshape(wout,size(xIMesh));
        %disp('size of WW='); size(WW)
        %figure
        surfc(xIMesh,yIMesh,WW)
        shading interp
        colorbar
        title({['\bfMode ',int2str(i)]; ['Eig=',num2str(EgN(i),'%9.2f')]})
    end

    % Plot modes shapes with no shading
    if( plot_flag == 1 || plot_flag == 3 )
        set(0,'CurrentFigure',h2);
        subplot(1,6,i)
        WW=reshape(wout,size(xIMesh));       
        M=mesh(xIMesh,yIMesh,WW);
        title({['\bfMode ',int2str(i)]; ['Eig=',num2str(EgN(i),'%9.2f')]})
    end
    
    % Plot animated shapes
    if( plot_flag == 2 )
        WW=reshape(wout,size(xIMesh));
        %disp('size of WW='); size(WW)
        figure
        h=surf(xIMesh,yIMesh,WW);
        shading interp
        axis([min(min(xIMesh)) max(max(xIMesh)) min(min(yIMesh)) max(max(yIMesh)) min(min(WW)) max(max(WW))]);
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
