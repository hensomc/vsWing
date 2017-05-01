%%  vMain
function [wingDATA,Lambda1] = main(varargin)
global Lambda1

% clear workspace, close all plots, clear command window
%clear % clear all variables
% close all
% clc
format short

fprintf('\n\n\n\n\n');
fprintf('============\n');
fprintf('%s\n', get_prog_vers);
fprintf('============\n');

degenGeomFileName='';

% Optional input argument
if length(varargin)==1
    jobID=varargin{1};
    close all
    clc
elseif length(varargin)>1
    jobID=0;
    degenGeomFileName=varargin{2};
else
    %clear
    close all
    clc
    % -----------------
    % Run a default job
    % -----------------
    %jobID=3;
    %jobID=4;
    %jobID=5; % Cantilevered plate uniform bending
    %jobID=8; % Cantilevered VAT2 laminate
    %jobID=20;
    %jobID=21; % Unswept cantilevered wing y-dir span, curv paths

    %jobID=97;  % Unswept cantilevered wing y-dir span (SPAR WEBS at midplane)
    %jobID=98;  % Unswept cantilevered wing-plate y-dir span (ribs & spars at midplane)
    %jobID=99;  % Unswept cantilevered wing y-dir span (SPAR CAPS at midplane)

    %jobID=100; % VSP swept wing
    %jobID=101 %  Cantilevered Trapezoidal Wing : Ritz to FEM comparison - Cantilevered with uniform pressure
    %jobID=102;  % Unswept cantilevered wing x-dir span
    %jobID=103;  % Unswept VSP wing y-dir span
    %jobID=104;  % Swept VSP wing y-dir span, spars
    %jobID=105;  % Unswept VSP wing y-dir,rect foil, spars
    %jobID=107;  % Unswept SOLID VSP wing y-dir,rect foil
    %jobID=108;  % Unswept SOLID Fill VSP Wing w/ 3SPAR CAPS y-dir span    

    %jobID=109;  % VSP F16 COARSE Wing   
    
    % =========
    % LIU CASES
    %jobID=106;  % Trapezoidal Plate Wing Liu, linear tapered thickness, isotropic
    %jobID=110;  % VSP LIU SOLID Wing 
    %jobID=111;   % VSP LIU Built-Up Wing, 4 Spars, 10 Ribs
    %jobID=112;   % VSP LIU BUILT-UP WING w/ SMEARED QUASI LAM, 4 SPARS, 10 RIBS
    %jobID=113;   % VSP LIU BUILT-UP WING w/ 4PLY SMEARED VLAM, 4 SPARS
    % =========
    
    % =======================
    % AIAA SciTech 2017 CASES
    %jobID=114;  % AIAA Uniform Trapezoidal Plate, linear tapered thickness, isotropic
    %jobID=115;  % AIAA Uniform Trapezoidal Core Filled Wing, Rect Foil
    %jobID=113;  % LIU BUILT-UP WING w/ 4PLY SMEARED VLAM, 4 SPARS
    % =======================

    jobID=96;   % Cantilevered Unswept Wing y-dir span


    % ===========
    % ESAVE CASES
    %jobID = 120;
    % ===========
    
    %jobID=200;  % Minimize max strain, vol constraint
    %jobID=201; % Maximize vib frequency, vol constraint
    %jobID=202; % Min weight, strain constraint
    %jobID=203; % Min weight, strain & buckling constraint

    % jobID=204; % Min weight, linear var_thk, strain constraint
    % jobID=205; % Min weight, quadratic var_thk, strain constraint

    % jobID=206; % Multi Obj: Min weight, Max Freq, strain constraint
    % jobID=207; % Multi Obj: Min weight, Max Freq, strain constraint, linear variable thickness
    % jobID=208; % Multi Obj: Min weight, Max Freq, strain constraint, quadratic variable thickness

    %jobID=209; % Multi Obj: Min weight, strain constraint, Trapezoidal Plate Wing Liu, linear tapered thickness, isotropic
end

% Get input data
[wingDATA,solDATA,lamDATA,ldsDATA,desDATA]=getInputs(jobID,degenGeomFileName);
%return;

% Run solution sequence
solStatus=solSequence(wingDATA,lamDATA,solDATA,ldsDATA,desDATA);

% Trial function to output all open plots
%publish('vlamWordOutput.m','doc')
%saveAllFigs2PDF();

% enables export of figures to powerpoint (use button icon on Fig dialog).
%pptfigure('all')


%% Run the solution sequence and return its status
function solStatus = solSequence(wingDATA,lamDATA,solDATA,ldsDATA,desDATA)

fprintf('Starting job, soltype = %i\n',solDATA.isoltype);
solStatus=0;

% Recover solution parameters
femSol=solDATA.femSol;
bcType=solDATA.bcType;

% Recover panel data
panXY = wingDATA.panXY;

% =============
% Analysis Only
% =============
if(solDATA.iDesOpt == 0)

    % Get Mesh
    %% Generate mesh data
    [ meshDATA, intDATA ] = meshGen( wingDATA, solDATA );
    [psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts,surfI,xyzI] = getIntData(intDATA);
    sMesh=psiIMesh;
    tMesh=etaIMesh;
    
    % Ritz Solution
    %[C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);
    [C,ritzWt,ritzEig,ritzDisp,ritzEps]=vLam(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA);
    
    % FEM Solution
    nasWt=0;nasEig=0;nasDisp=0;nasEps=0;XYData1=0;
    if(femSol==1)
        [nasWt,nasEig,nasDisp,nasEps,nasStatus,XYData1]=wingFEM(wingDATA,solDATA,lamDATA,ldsDATA);
    end
    
    % Post Processing
    getSolErr(panXY,solDATA,lamDATA,C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh,nasWt,nasEig,nasDisp,nasEps,XYData1);
    
elseif(solDATA.iDesOpt == 1)
    % ===================
    % Design optimization
    % ===================
    % Adjust solution parameters for speed
    solDATA.pDeg=4;       
    %solDATA.Mg=30;
    solDATA.Mg=10;
    solDATA.femSol=0;       % 1=create a Nastran FE model for solution
    solDATA.plot_flag=0;    % 0 = none, 1=modes, 2=animate modes, 3=details
    solDATA.iThkFun=1;      % 1=use polynomial thickness function

    Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA);

    % ===================
    % === PARAM STUDY ===
    % ===================    
elseif( solDATA.iDesOpt == 2 )
    
    % Set parameter study type
    %paramType = 2;  % 1 = Polynomial Order, 2= Integ Pts Mg, 3 = Plate Skew, 4 = Layup Angle <T0,T1>
    paramType = solDATA.paramType;

    % Defaults: outer paramter loop
    pval2I=0; pval2F=0; pval2D=1;
    
    switch paramType
        case 1  % polynomial order
            pvalI=solDATA.pDeg;
            pvalF=pvalI; pvalD=1;
            pvalI=4; pvalF=12; pvalD=1;

            
        case 2  % integration order
            pvalI=solDATA.Mg;
            pvalF=pvalI;
            pvalI=10; pvalF=40; pvalD=5;

            
        case 3  % panel skew angle
            pvalI=0;
            pvalF=pvalI;
            pvalI=0; pvalF=60; pvalD=15;

            
        case 4  % layup angle <T0|T1>
            % <T0=45|T1=-90:90> : Houmat laminate
            pvalI=45; pvalF=-90; pvalD=-15;
            
            % <T1=0:90) : Gurdal Case I laminate
            pvalI=0; pvalF=90; pvalD=15;
            
            % <T0=0:90> outer paramter loop
            pval2I=0; pval2F=90; pval2D=15;
            
            
            
            % <Phi=90|T1=0:90) : Gurdal Case II laminate
            % <T1=0:90) : Gurdal Case II laminate
            %pvalI=0; pvalF=90; pvalD=15;
            
            % <T0=90:0> outer paramter loop
            %pval2I=90; pval2F=0; pval2D=-15;
    end
    paramStudy( pvalI, pvalF, pvalD, pval2I, pval2F, pval2D, paramType, panXY,lamDATA,solDATA,ldsDATA);
      
else
    disp('ERROR: iDesOpt case not supported: '); solDATA.iDesOpt
    return;  
end



function paramStudy( pvalI, pvalF, pvalD, pval2I, pval2F, pval2D, paramType, panXY,lamDATA,solDATA,ldsDATA );

%% Recover solution parameters
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,~,numMode,ssCalc,paramType]=getSolData(solDATA);

if(plot_flag == 1)
    pvalI
    pvalF
    pvalD
end
    
    panXY_orig=panXY;

% Initialize graphing variables
xval=[];
yval=[];
xvalFEM=[];
yvalFEM=[];


%% Loop over parameter values

% Outer parameter loop
j2=0;
for pval2=pval2I:pval2D:pval2F
    pval2
    j2=j2+1;
    
    xval=[]; yval=[]; yvalFEM=[]; ratio=[]; ratio1=[];
    %clear freqParam_ritz1 freqParam_nas;
    
    % Inner parameter loop
    i=0;
    for pval=pvalI:pvalD:pvalF
        i=i+1;
        
        pval
        
        % Set solution parameters
        switch paramType
            case 1  % polynomial order
                solDATA.pDeg = pval;
                if(solDATA.pDeg > solDATA.Mg)
                    solDATA.Mg=solDATA.pDeg;
                end
                
            case 2  % integration order
                solDATA.Mg = pval;
                
            case 3  % panel skew angle
                % Modify position of panel corners 3 & 4: initially rectangular
                b=panXY_orig(3,2)-panXY_orig(2,2);
                x=b*sind(pval);
                y=b*cosd(pval);
                panXY(3,1)=panXY_orig(3,1)+x; panXY(3,2)=y;
                panXY(4,1)=panXY_orig(4,1)+x; panXY(4,2)=y;
                
            case 4  % layup angle <T0|T1>
                
                T0=lamDATA.theta0(1);  % default is no change
                if(pval2F ~= pval2I)
                    T0=pval2;
                end
                
                T1=pval;
                %setVarLamTheta( lamDATA, T0, T1 );
                nply=length(lamDATA.theta0);
                for iply=1:nply/2
                    if(mod(iply,2) == 1)
                        lamDATA.theta0(iply)=T0;
                        lamDATA.theta1(iply)=T1;
                    else
                        lamDATA.theta0(iply)=-T0;
                        lamDATA.theta1(iply)=-T1;
                    end
                end
                j=1;
                for iply=nply:-1:nply/2+1
                    lamDATA.theta0(iply)=lamDATA.theta0(j);
                    lamDATA.theta1(iply)=lamDATA.theta1(j);
                    j=j+1;
                end
                
                % Compute effective laminate properties for normalization
                [eqx, eqy, eqxn, eqyn] = get_equiv_lam_props(panXY, lamDATA, solDATA.smearKM);
        end

               
        % Ritz Solution
        %[~,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
        
        % Generate mesh data
        [ meshDATA, intDATA ] = meshGen( wingDATA, solDATA );
        [psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
        [psiIMesh,etaIMesh,psiI,etaI,gl_wts,surfI,xyzI] = getIntData(intDATA);
        sMesh=psiIMesh;
        tMesh=etaIMesh;
        
        % Ritz Solution
        %[C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);
        [C,ritzWt,ritzEig,ritzDisp,ritzEps]=vLam(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA);

        % FEM Solution
        if(femSol==1)
            nasWtVals=0;nasEigVals=0; nasDisp=0; XYData1=0;
            [nasWtVals,nasEigVals,nasDisp,nasEps,nasStatus,XYData1]=wingFEM(wingDATA,solDATA,lamDATA,ldsDATA);
        end
        
        % Gather results for plot
        if(isoltype == 1 || isoltype == 2)  % eigenvalue (vib or buckling)
            if(length(ritzEig) <= numMode)
                xval(i,1:length(ritzEig))=pval; yval(i,1:length(ritzEig))=ritzEig(1:length(ritzEig));
                %xval(i,1:length(ritzEig))=eqxn(1);  % Temporary hardwired values
                if(femSol==1)
                    xvalFEM=xval;
                    yvalFEM(i,1:length(ritzEig))=nasEigVals(1:length(ritzEig));
                end
            else
                xval(i,1:numMode) = pval; yval(i,1:numMode)=ritzEig(1:numMode);
                %xval(i,1:length(ritzEig))=eqxn(1); % Temporary hardwired values
                if(femSol==1)
                    xvalFEM=xval; yvalFEM(i,1:numMode)=nasEigVals(1:numMode);
                    ratio(i,1:numMode)=yval(i,1:numMode)./yvalFEM(i,1:numMode);
                end
            end
            
        elseif(isoltype == 3)  % deflection
            xval(i,1)=pval;
            if( max(any(abs(ldsDATA.ptLds))) > 0 || any(abs(ldsDATA.p0)) > 0 )
                yval(i,1)=max(max(abs(ritzDisp(:,3))));
            else
                %yval(i,1)=max(max(abs(Vout)));
                yval(i,1)=max(max(abs(ritzDisp(:,2))));
            end
            if(femSol==1)
                xvalFEM=xval; yvalFEM(i,1)=max(max(abs(nasDisp(:,2))));
            end
        end
        
        % Fundamental frequency parameter
%         if(isoltype == 1 || isoltype == 2)
%             [freqParam_ritz1(i)]=getFreqParam(ritzEig(1),panXY,lamDATA,1);
%             [freqParam_ritz2(i)]=getFreqParam(ritzEig(1),panXY,lamDATA,2);
%             if(femSol==1)
%                 [freqParam_nas1(i)] =getFreqParam(nasEigVals(1),panXY,lamDATA,1);
%                 [freqParam_nas2(i)] =getFreqParam(nasEigVals(1),panXY,lamDATA,2);
%             end
%         end
    end
    
    % Post process
    paramOut(panXY, solDATA, lamDATA, paramType, pvalI, pvalF, pval2I, pval2F, pval2, j2, xval, yval, xvalFEM, yvalFEM, ratio);
    

end
% end of outer paramter loop
if(pval2F~=pval2I)
    hold off
end

function paramOut(panXY, solDATA, lamDATA, paramType, pvalI, pvalF, pval2I, pval2F, pval2, j2, xval, yval, xvalFEM, yvalFEM, ratio)
%% Tabulate
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,~,numMode,ssCalc,paramType]=getSolData(solDATA);

persistent ptxt;

% Tabluate parametric results
if( isoltype == 1 || isoltype == 2)
    %tritz=table(xval(:,1), yval(:,1), yval(:,2), yval(:,3), yval(:,4), yval(:,5));
    %tritz=table(xval(:,1), yval(:,1), yval(:,2), yval(:,3), yval(:,4), yval(:,5));
    tritz=table(xval(:,1),yval(:,1:numMode));
    if( femSol==1)
        tnas=table(xval(:,1), yvalFEM(:,1), yvalFEM(:,2), yvalFEM(:,3), yvalFEM(:,4), yvalFEM(:,5));

        t1=table(xval(:,1), yval(:,1), yvalFEM(:,1), ratio(:,1));
        tr=table(xval(:,1), ratio(:,1), ratio(:,2), ratio(:,3), ratio(:,4), ratio(:,5));
    else
        %t1=table(xval(:,1), yval(:,1), yval(:,2), yval(:,3), yval(:,4), yval(:,5));
        t1=table(xval(:,1),yval(:,1:numMode));
    end
    disp('Summary of Ritz Eigenvalues');
    disp(tritz);
    
    if(femSol==1)
        disp('Summary of Nas Eigenvalues');
        disp(tnas);
    end
    
    disp('Summary of 1st Eigenvalue');
    disp(t1)
    
    if(femSol==1)
        disp('Summary of 1st 5 Eigenvalue Ratios');
        disp(tr)
    end
    
    % Tabluate frequency parameters
%     if( isoltype == 1 || isoltype == 2)
%         if(femSol==1)
%             ratio1=freqParam_ritz1(:)./freqParam_nas1(:);
%             ratio2=freqParam_ritz2(:)./freqParam_nas2(:);
%             t2=table(xval(:,1), freqParam_ritz1(:), freqParam_nas1(:), ratio1, freqParam_ritz2(:), freqParam_nas2(:), ratio2);
%         else
%             t2=table(xval(:,1), freqParam_ritz1(:), freqParam_ritz2(:));
%         end
%         disp('Summary of 1st Mode Non-Dimensional Frequency Parameter');
%         disp(t2)
%     end
end

%% Parameter Plots
if( pvalF ~= pvalI )
    
    % Color matrix for multiple curves (10 Max)
    if(pval2F ~= pval2I)
        scm = ['r-';      % red solid
            'g:';      % green dotted
            'b-';      % blue solid
            'm:';      % magenta dotted
            'c-';      % cyan solid
            'y:';      % yellow dotted
            'r:';      % red dotted
            'g-';      % greeb solid
            'b:';      % blue dotted
            'y-'];     % yellow solid
    end
    
    % Get plot titles
    [panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);
    panTitle_lamTitle = strcat(panTitle,',  ',lamTitle);
    
    % Frequency convergence
    if(isoltype == 1 || isoltype ==2)
        if(pval2F~=pval2I && pval2==pval2I)
            figure;
        elseif(pval2F==pval2I)
            figure;
        end
        
        % Plot ritz data
        if(pval2F~=pval2I)
            plot(xval,yval,scm(j2,:),'linewidth',2);
        else
            plot(xval,yval,'-o','linewidth',2);
        end
        
        % Plot fem data
        if(femSol==1)
            hold on
            plot(xvalFEM,yvalFEM,':s','linewidth',1);
        end
        
        % Set title
        if(isoltype == 1)
            title( {'\bfVibration Eigenvalue Convergence';panTitle_lamTitle;''},'fontsize',12 )
        else
            title( {'\bfBuckling Eigenvalue Convergence';panTitle_lamTitle;''},'fontsize',12 )
        end
        
        ylabel('\bfEigenvalue','fontsize',12);
        % Build legend
        if(pval2F==pval2I)
            for i=1:numMode
                if( femSol==1 )
                    txt(i)={['Ritz Mode ',num2str(i)]};
                else
                    txt(i)={['Mode ',num2str(i)]};
                end
            end
        else
            %txt(j2)={['pval2: ', num2str(pval2)]};
            txt=strcat('T0: ', num2str(pval2));
            %ptxt{j2}={['pval2: ', num2str(pval2)]};
            ptxt{j2}=txt;
        end
        
        % Legend items for FEM results
        if( femSol==1 )
            j=0;
            for i=numMode+1:2*numMode
                j=j+1;
                txt(i)= {['FEM Mode ',num2str(j)]};
            end
        end
        
        % Draw legend
        if(pval2F==pval2I)
            legend(txt,'location','northeast','FontSize',8);
            axis([pvalI,pvalF,0, max(max(yval))])
        else
            if(pval2==pval2F)
                %legend(txt,'location','best','FontSize',8);
                legend(ptxt,'location','best','FontSize',8);
            end
        end
        
    % Displacement convergence
    elseif(isoltype == 3)
        if(pval2F~=pval2I && pval2==pval2I)
            figure;
        end
        plot(xval,yval,'r-.o');
        if(femSol==1)
            hold on
            plot(xvalFEM,yvalFEM,'b:s');
        end

        title( {'\bfMax Deflection Convergence'},'fontsize',14 )
        ylabel('\bfDeflection (in.)','fontsize',12);
        txt(1)={['Ritz Max Displ']};
        if(femSol==1)
            txt(2)={['FEM Max Displ']};
        end
        legend(txt,'location','northeast','FontSize',8,'location','best');
        axis([pvalI,pvalF,0, max(max(yval))])
    end
    
    % Xlabel dependent on type        
    switch paramType
        case 1
            xlabel('\bfRitz Polynomial Order (N_{p})','fontsize',12)
        case 2
            xlabel('\bfIntegration Order (M_{g}xN_{g})','fontsize',12)
        case 3
            xlabel('\bfPlate Skew Angle \phi_{skew} (Deg)','fontsize',12)
            title({'\bfVibration Eigenvalues for Skew Panels';panTitle_lamTitle;''},'fontsize',14 )
        case 4
            if(isoltype==1)
                title( {'\bfVibration Eigenvalues for VAT Laminates',panTitle},'fontsize',14 )
            elseif(isoltype==2)
                title( {'\bfBuckling Factors for VAT Laminates',panTitle},'fontsize',14 )
            end
            xlabel('\bfLayup Angle T1 (Deg)','fontsize',12)
            %xlabel('\bfNormalized Stiffness E_x^{eq}/E_1','fontsize',12)
            
    end
end

if(pval2F ~= pval2I)
    hold on
    axis auto
end
