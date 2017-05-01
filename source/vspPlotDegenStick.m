function [degenGeom,panXY] = plotDegenStick( dgfile, plot_flag )

% Load VSP output degenerate geometry.
run(dgfile);

% Initialize flags to prevent undefined access.
plotLE = false; plotTE = false; plotCGSolid = false; plotCGShell = false;
plotInt = false; plotThick = false;

if plot_flag>0
    plotLE = true;
    plotTE = true;
    plotInt = true;
    plotCGSolid = true;
    plotCGShell = true;
    plotThick = true;
end

intpt = 0.25;

ngeom = length(degenGeom);
disp({'ngeom = ' ngeom});

panXY=[];

%figure(1)
if plot_flag>0
    figure
    clf
    hold on
end

for i=1:ngeom
  disp(['Component ' num2str(i)  ' Name: ' degenGeom(i).name  ' Type: ' degenGeom(i).type]);

  for j=1:length(degenGeom(i).stick)
%     disp(['length(degenGeom(i).stick = ' int2str(length(degenGeom(i).stick))]);
%     disp(['length(degenGeom(i).stick(j).lex = ' int2str(length(degenGeom(i).stick(j).lex))]);
%     disp(['length(degenGeom(i).stick(j).ley = ' int2str(length(degenGeom(i).stick(j).ley))]);
%     disp(['length(degenGeom(i).stick(j).lez = ' int2str(length(degenGeom(i).stick(j).lez))]);

    xle = degenGeom(i).stick(j).lex;
    yle = degenGeom(i).stick(j).ley;
    zle = degenGeom(i).stick(j).lez;

    xte = degenGeom(i).stick(j).tex;
    yte = degenGeom(i).stick(j).tey;
    zte = degenGeom(i).stick(j).tez;

    xint = xle + intpt*(xte-xle);
    yint = yle + intpt*(yte-yle);
    zint = zle + intpt*(zte-zle);

    if plot_flag>0
        if(plotLE)
            plot3(xle, yle, zle, 'kx-')
        end
        
        if(plotTE)
            plot3(xte, yte, zte, 'bo-')
        end
        
        if(plotInt)
            plot3(xint, yint, zint, 'k-.')
        end
        
        if(plotCGSolid)
            plot3(degenGeom(i).stick(j).cgSolidx,degenGeom(i).stick(j).cgSolidy,degenGeom(i).stick(j).cgSolidz,'r--');
        end
        
        if(plotCGShell)
            plot3(degenGeom(i).stick(j).cgShellx,degenGeom(i).stick(j).cgShelly,degenGeom(i).stick(j).cgShellz,'y--');
        end
        
        if(plotThick)
            xt = xle + degenGeom(i).stick(j).tLoc.*(xte-xle);
            yt = yle + degenGeom(i).stick(j).tLoc.*(yte-yle);
            zt = zle + degenGeom(i).stick(j).tLoc.*(zte-zle);
            
            % Max thickness line
            plot3(xt,yt,zt,'k--');
            
            % Needs normal from plate.
            dx = 0.5 * degenGeom(i).stick(j).toc .* degenGeom(i).stick(j).chord .* degenGeom(i).plate(j).nx;
            dy = 0.5 * degenGeom(i).stick(j).toc .* degenGeom(i).stick(j).chord .* degenGeom(i).plate(j).ny;
            dz = 0.5 * degenGeom(i).stick(j).toc .* degenGeom(i).stick(j).chord .* degenGeom(i).plate(j).nz;
            
            % Points above and below max thickness.
            plot3(xt+dx,yt+dy,zt+dz,'kx',xt-dx,yt-dy,zt-dz,'kx');
            
        end
        hold off
        
        axis equal
        axis off
    end
    
    
    % Assign plate coordinates (i=1, RHS, assumes only 1 component)
    if(i == 1)
        numseg = (length(degenGeom(i).stick(j).lex)-1)/5;
        indx=1;
        for k=1: numseg
            segXY = [xle(indx) yle(indx);
                xte(indx) yte(indx);
                xte(indx+5) yte(indx+5);
                xle(indx+5) yle(indx+5)] * 12;       
    %             xte(length(xte)) yte(length(yte));
    %             xle(length(xle)) yle(length(yle))] * 12  
            indx=indx+5;
            panXY{k}=segXY;
        end

        
        %transform wing - temporary fix to address vs_lam solve issue
%         panXY = [panXY(1,1) panXY(1,2);
%             panXY(4,2) panXY(4,1);
%             panXY(3,2) panXY(3,1);
%             panXY(2,2) panXY(2,1)];
    end
             
  end
end