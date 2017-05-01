function [ output_args ] = getSolErr( panXY,solDATA,lamDATA,C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh,nasWt,nasEig,nasDisp,nasEps,XYData1 )

%getSolErr Compares Ritz solution with FEM solution
%   Detailed explanation goes here

% Recover solution parameters
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);


%% Compute and Plot Modes/Deformed Shape
if(isoltype == 1 || isoltype == 2)
    plotModes(panXY,ipoltype,pDeg,bcMethod,sMesh,tMesh,C,ritzEig,plot_flag)
elseif(isoltype == 3)
    plotDefShape(panXY,solDATA,lamDATA,sMesh,tMesh,C,plot_flag);
end

%% Weight Summaries
disp('======================');
disp('Ritz Weight Summary');
disp('======================');
ritzWt

disp('======================');
disp('Nastran Weight Summary');
disp('======================');
nasWt

disp('======================');
disp('panelWeightVthk Weight');
disp('======================');
panWtVthk = panelWeightVthk(solDATA,lamDATA,panXY )

%% Compute error in solution
if(isoltype == 1 || isoltype == 2)
    if(femSol == 1)
        eigError=getEigErrFEM(panXY,solDATA,lamDATA,nasEig,ritzEig,plot_flag,numMode);
    else
        if(plot_flag>0)
            disp('Summary of Ritz Eigenvalues');ritzEig(1:numMode)
        end
    end
    
    % === Static Solution: Displacement and strain error ===
elseif(isoltype == 3 && femSol == 1)
    
    getDispErrFEM(panXY,XYData1,solDATA,lamDATA,nasDisp,ritzDisp,plot_flag);
    if(ssCalc == 1 )
        getEpsErrFEM(panXY,XYData1,solDATA,lamDATA,nasEps,ritzEps,plot_flag);
    end
end
