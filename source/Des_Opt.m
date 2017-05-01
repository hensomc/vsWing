function [XOPT,his]= Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% Design Optimization - Main driver function
disp('DES_OPT called');

global h_histFig;
h_histFig=0;

plot_hist_flag=1;
numX = length(desDATA.X0);

%desResults=cell(length(desDATA.Xstart));
desResults=[];

% Generate mesh data
[ meshDATA, intDATA ] = meshGen( wingDATA, solDATA );

% Loop over start points
disp('length(desDATA.Xstart='); length(desDATA.Xstart)
for i=1:length(desDATA.Xstart)
    desDATA.X0 = desDATA.Xstart{i};

    switch desDATA.drespType
        case {'WEIGHT'}
            [XOPT,his]= Weight_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA,desDATA);
        case {'EIGN'}
            [XOPT,his]= Vib_Freq_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
        case {'BUCKL'}
            [XOPT,his]= Buckl_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
        case {'STRAIN'}
            [XOPT,his]= Strain_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
        case {'WEIGHT_EIGN'}
            [XOPT,his,hisF]= Weight_Eign_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
            numX=length(XOPT);
            [his hisF]
            his=[his hisF(:,numX+1)]
        otherwise
            return;
    end
        
    % Disable plotting to dynamic history plot
    if gcf==h_histFig
        hold off;
        %disp('hold off of dyn hist plot');
    end
    
    %desResults{i} = {XOPT his(numX+1,end) his(numX+2,end)};
    %desResults = [desResults; XOPT(1:end) his(end,numX+1) his(end,numX+2)];
    desResults = [desResults; XOPT(1:end) his(end,numX+1:end)];

    % Plot full history results - redundant with dynamic history plot
    plotDes(XOPT,his,wingDATA,lamDATA,solDATA,ldsDATA,desDATA,plot_hist_flag,i);
    
end
% celldisp(desResults)
% t1=cell2table(desResults);
% disp(t1)
%format long
format short
desResults

% Find optimum design values for all start points
disp('Iteration with optimum');
I=find(max(max(desResults(:,numX+1))))  % max good only for max Pcrit?
XOPT=desResults(I,1:numX)

disp('Pausing to allow visual check, hit any key to continue...')
%pause

% Check optimal solution (Ritz and FEM)
solDATA.pDeg=8;
%solDATA.Mg=solDATA.pDeg;
solDATA.Mg=30;
sz=size(wingDATA.spar.XY);
nSpar=sz(2);
solDATA.Mg = check_Mg(nSpar,solDATA.Mg);
solDATA.femSol=1;
solDATA.plot_flag=1;

% ... Regenerate mesh
[ meshDATA, intDATA ] = meshGen( wingDATA, solDATA );
sMesh=intDATA.psiIMesh;
tMesh=intDATA.etaIMesh;

% Update laminate design values
nvar=0;
% if(lamDATA.tDeg>0)
%     lamDATA.thk_coeff = XOPT;
%     % Draw surface thickness constraint
%     [g,heq]=vlamThkCon(XOPT,wingDATA,lamDATA,solDATA,meshDATA,intDATA);
%     %pause
% else
    switch desDATA.drespType
        case {'WEIGHT' 'EIGN'}
            
            % loop over dvType
            for i=1:length(desDATA.dvType)
                
                % =====================
                % Thickness design Vars
                % =====================
                if desDATA.dvType(i)==0
                    for j=1:length(lamDATA.thk)
                        nvar=nvar+1;
                        lamDATA.thk(j)=XOPT(nvar);
                    end
                end
                
                % =======================
                % Fiber angle design Vars
                % =======================
                % Fiber angles (i.e. T0, T1), link +45,-45,90 plies
                if desDATA.dvType(i)==1
                    ntheta=2; % hardwired for linear theta variation
                    lamDATA.theta0(1)=XOPT(nvar+1);
                    lamDATA.theta1(1)=XOPT(nvar+2);
                    nvar=nvar+ntheta;
                    
                    lamDATA.theta0(2)=lamDATA.theta0(1)+45;
                    lamDATA.theta0(3)=lamDATA.theta0(1)-45;
                    lamDATA.theta0(4)=lamDATA.theta0(1)+90;

                    lamDATA.theta1(2)=lamDATA.theta1(1)+45;
                    lamDATA.theta1(3)=lamDATA.theta1(1)-45;
                    lamDATA.theta1(4)=lamDATA.theta1(1)+90;
                end
                
                % =================================
                % Polynomial Thickness Coefficients
                % =================================
                if desDATA.dvType(i)==2
                    for j=1:length(lamDATA.thk_coeff)
                        nvar=nvar+1;
                        lamDATA.thk_coeff(j)=XOPT(nvar);
                    end
                end
            end
            
        case {'BUCKL' 'STRAIN'}
            % Fiber angles only
            %setVarLamTheta( lamDATA, round(XOPT(1)), round(XOPT(2)) );
            
            nply=length(lamDATA.theta0);
            for i=1:nply/2
                if(mod(i,2) == 1)
                    lamDATA.theta0(i)=round(XOPT(1));
                    lamDATA.theta1(i)=round(XOPT(2));
                else
                    lamDATA.theta0(i)=-round(XOPT(1));
                    lamDATA.theta1(i)=-round(XOPT(2));
                end
            end
            j=1;
            for i=nply:-1:nply/2+1
                lamDATA.theta0(i)=lamDATA.theta0(j);
                lamDATA.theta1(i)=lamDATA.theta1(j);
                j=j+1;
            end
            
        otherwise
            return;
    end
% end

% Compare Ritz optimized solution with Nastran model
%[C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);
[C,ritzWt,ritzEig,ritzDisp,ritzEps]=vLam(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA);
disp('finished vLam call');
[nasWt,nasEig,nasDisp,nasEps,nasStatus,XYData1]=wingFEM(wingDATA,solDATA,lamDATA,ldsDATA);
disp('finished wingFEM call');
%getSolErr(wingDATA.panXY,solDATA,lamDATA,C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh,nasWt,nasEig,nasDisp,nasEps,XYData1);
getSolErr(wingDATA.panXY,solDATA,lamDATA,C,ritzWt,ritzEig,ritzDisp,ritzEps,sMesh,tMesh,nasWt,nasEig,nasDisp,nasEps,XYData1);
disp('finished getSolErr call');

        
function [XOPT,his]= Weight_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA,desDATA)
%% Weight_Des_Opt.m : Min weight design optimization w/ buckling constraint
disp('Min weight design optimization w/ strain and/or buckling constraint');

global his;
his=[];

global Lambda1;

global hisG;
hisG=[];

global PlotThkSurf;
numIter=100;  % max optim iterations
numX=10;
nfval=numX*numIter;
PlotThkSurf=zeros(nfval);
global iPlotThkSurf;
iPlotThkSurf=0;

% Handles to objective and constraint functions
OBJ=@Weight_Des_Obj;
CON=@genCon;

% Thickness inequality constraints
A=[];b=[];
if( lamDATA.tDeg > 0 )
%     [A]=getInequalLayerCon(wingDATA,lamDATA,solDATA,meshDATA,intDATA);
% %      tmax(size(A))= 0.150;
%       tmin(size(A))= 0.020;
%       tmin=tmin(:);
% 
% %   Min thk constraint
%     A=-A';
%     b=-tmin;
%     
% %   Max & Min thk constraint
% %     A=[A';-A'];
% %     b=[tmax;-tmin];
end

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
%disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
%OP=[];
%OP=optimset('disp','iter')  % Display optimization/iteration details
%XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA)
%XOPT=fmincon(OBJ,X0,A,b,[],[],XL,XU,CON,OP,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)

% Set optimizer options
options=optimoptions('fmincon');
%options.TolFun=1.0e-3;
options.MaxIter=100;
% options.TolX=1.0e-4;
%options.TolCon=1.0e-3;
options.Display='iter';
options.UseParallel=true;
% options.PlotFcns={@optimplotfval};
% options.OutputFcn=@outfun;

% output function
%options = optimoptions(@fmincon,'OutputFcn',@outfun,'Display','iter','Algorithm','active-set');

% plot functions
%options=optimset('PlotFcns',@optimplotfval);
%options = optimset('PlotFcns',{@optimplotx,@optimplotfval,@optimplotfirstorderopt});

disp('Optimizer options='); options

%XOPT=fmincon(OBJ,X0,A,b,[],[],XL,XU,CON,options,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
[XOPT,FVAL,ExitFlag,Output,Lambda1]=fmincon(OBJ,X0,A,b,[],[],XL,XU,CON,options,wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA,desDATA)

% Check constraint
if length(A)>0
    tt=A*XOPT.';
    tt_min=min(tt)
end

% Place design and constraint history into one matrix
his= [his hisG];



function [XOPT,his]= Vib_Freq_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_OPT.m : Max frequency design optimization w weight/vol constraint
disp('Optimize plate for max frequency subject to volume constraint');

global his
his=[];

global hisG;
hisG=[];

% Handles to objective and constrint functions
OBJ=@vlam_Freq_Obj;
CON=@vlam_Volume_Con;
CON=[];
% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0
%his = X0;  % Initial values X0 not captured by "his" - this fails out

% Inequality layer-thickness constraints
if( lamDATA.tDeg > 0 )
    [A]=getInequalLayerCon(wingDATA,lamDATA,solDATA)
    %pause
%     A = [A;-A];
     tmax(size(A))= 0.150;
     tmin(size(A))= 0.001;
%     b = [tmax;-tmin];

%   Min thk constraint
    A=-A';
    b=-tmin;
    
%   Max & Min thk constraint
%     A=[A';-A'];
%     b=[tmax;-tmin];
end

% Optimize
OP=[];
OP=optimset('disp','iter');  % Display optimization/iteration details
XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
%XOPT=fmincon(OBJ,X0,A,b,[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA)

% Place design and constraint history into one matrix
his= [his hisG];

function [XOPT,his]= Buckl_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_OPT.m : Max buckling eigenvalue design optimization w weight/vol constraint
disp('Optimize plate for max buckling eigenvalue subject to volume constraint');

global his
his=[];

% Handles to objective and constrint functions
OBJ=@Buckl_Obj;
CON=@fibK_Con;

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
OP=optimset('disp','iter')  % Display optimization/iteration details
XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)


function [XOPT,his]= Strain_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_OPT.m : Minimize Max strain design optimization w weight/vol constraint
disp('Minimize Max strain design optimization w weight/vol constraint');

global his
his=[];

% Handles to objective and constrint functions
OBJ=@Strain_Obj;
%CON=@fibK_Con;
CON=@vlam_Volume_Con;

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
OP=optimset('disp','iter')  % Display optimization/iteration details
XOPT=fmincon(OBJ,X0,[],[],[],[],XL,XU,CON,OP,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)


%function [f,hx,xC]=Weight_Des_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
function [f]=Weight_Des_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA,desDATA)
%% Weight_Des_Obj.m - weight design objective function g

global his;

%if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness
    %     t45=max(XDES(2), XDES(3));
    %     lamDATA.thk(1)=XDES(1);
    %     lamDATA.thk(2)=t45;
    %     lamDATA.thk(3)=t45;
    %     lamDATA.thk(4)=XDES(4);
    nvar=0;
    for i=1:length(desDATA.dvType)
        if desDATA.dvType(i)==0
            for j=1:length(lamDATA.thk)
                nvar=nvar+1;
                lamDATA.thk(j)=XDES(nvar);
            end
        end
       
        if desDATA.dvType(i)==1
            ntheta=2; % hardwired for linear theta variation
            lamDATA.theta0(1)=XDES(nvar+1);
            lamDATA.theta1(1)=XDES(nvar+2);
            nvar=nvar+ntheta;
            
            lamDATA.theta0(2)=lamDATA.theta0(1)+45;
            lamDATA.theta0(3)=lamDATA.theta0(1)-45;
            lamDATA.theta0(4)=lamDATA.theta0(1)+90;
            
            lamDATA.theta1(2)=lamDATA.theta1(1)+45;
            lamDATA.theta1(3)=lamDATA.theta1(1)-45;
            lamDATA.theta1(4)=lamDATA.theta1(1)+90;
        end
        
        %lamDATA.thk=XDES(1:2); % constant thickness
        %     lamDATA.theta0(1)=XDES(3);  % fiber angles
        %     lamDATA.theta1(1)=XDES(4);
        %     lamDATA.theta0(2)=-XDES(3);
        %     lamDATA.theta1(2)=-XDES(4);
        
        % =================================
        % Polynomial Thickness Coefficients
        % =================================
        if desDATA.dvType(i)==2
            for i=1:length(lamDATA.thk_coeff)
                nvar=nvar+1;
                lamDATA.thk_coeff(i)=XDES(nvar);
            end
        end
    end
% else
%     lamDATA.thk_coeff=XDES; % variable thickness
% end

f = panelWeightVthk(solDATA,lamDATA,wingDATA.panXY );

% Get current constraint value
%[g,heq]=buckl_Con(XDES,panXY,lamDATA,solDATA,ldsDATA);
% [C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(panXY,lamDATA,solDATA,ldsDATA);
% g=-EgN(1)+1;

% Add XDES, Weight and eigenvalue to design history
%his=[his;XDES f g];
his=[his;XDES f ];


function [g,heq]=buckl_Con(XDES,wingDATA,lamDATA,solDATA,ldsDATA)
%% buckl_Con.m : buckling constraint function

global his;

% if(lamDATA.tDeg==0)
%     %lamDATA.thk=XDES; % constant thickness
%     lamDATA.thk=XDES(1:2); % constant thickness
% %     lamDATA.theta0(1)=XDES(3);  % fiber angles
% %     lamDATA.theta1(1)=XDES(4);
% %     lamDATA.theta0(2)=-XDES(3);
% %     lamDATA.theta1(2)=-XDES(4);
% else
%     lamDATA.thk_coeff=XDES; % variable thickness
% end

% Buckling eigenvalue solution
solDATA.isoltype=3;
[C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);
g=-EgN(1)+1;

heq=0;


function [g,heq]=genCon(XDES,wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA,desDATA)
%% genCon.m : general constraint function

global his;

global hisG;

global PlotThkSurf;
global iPlotThkSurf;

% if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness
    nvar=0;
    for i=1:length(desDATA.dvType)
        if desDATA.dvType(i)==0
            for j=1:length(lamDATA.thk)
                nvar=nvar+1;
                lamDATA.thk(j)=XDES(nvar);
            end
        end        
        
        if desDATA.dvType(i)==1
            ntheta=2; % hardwired for linear theta variation
            lamDATA.theta0(1)=XDES(nvar+1);
            lamDATA.theta1(1)=XDES(nvar+2);
            nvar=nvar+ntheta;
            
            lamDATA.theta0(2)=XDES(1)+45;
            lamDATA.theta0(3)=XDES(1)-45;
            lamDATA.theta0(4)=XDES(1)+90;
            
            lamDATA.theta1(2)=XDES(2)+45;
            lamDATA.theta1(3)=XDES(2)-45;
            lamDATA.theta1(4)=XDES(2)+90;
        end
        
        %lamDATA.thk=XDES(1:2); % constant thickness
        %     lamDATA.theta0(1)=XDES(3);  % fiber angles
        %     lamDATA.theta1(1)=XDES(4);
        %     lamDATA.theta0(2)=-XDES(3);
        %     lamDATA.theta1(2)=-XDES(4);
        
        % =================================
        % Polynomial Thickness Coefficients
        % =================================
        if desDATA.dvType(i)==2
            for i=1:length(lamDATA.thk_coeff)
                nvar=nvar+1;
                lamDATA.thk_coeff(i)=XDES(nvar);
            end
        end
    end
% else
%     lamDATA.thk_coeff=XDES; % variable thickness
% end

% Static solution
g_strain=[];
solDATA.isoltype=3;
%[C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);
[C,ritzWt,EgN,ritzDisp,ritzEps]=vLam(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA);
% %ritzEps
[g_strain,heq_strain]=strain_Con(ritzEps,wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA);
%g_strain
%pause

% FEM solution
%[Wt,Eig,Disp,Eps,nasStatus,XYData1]=wingFEM(wingDATA,solDATA,lamDATA,ldsDATA);
%Eps(:,4)=1.0e-6;
%[g_strain,heq_strain]=strain_Con(Eps,wingDATA,lamDATA,solDATA,ldsDATA);
%g_strain
%pause;
% Buckling eigenvalue solution
g_buckl=[];
% solDATA.isoltype=2;
% [C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);
% g_buckl=-EgN(1)+1;

% Thickness constraints - Only need this constraint applied once (done)
g_thk=[];heq_thk=[];
if(lamDATA.tDeg>0)
    [g_thk,heq_thk]=vlamThkCon(XDES,wingDATA,lamDATA,solDATA,meshDATA,intDATA);
end
%disp('size of g_thk='); size(g_thk)
g=[g_buckl; g_strain; g_thk];
%g=[g_buckl; g_strain];
%disp('size of g='); size(g)

% Inequality constraints
heq=0;
% ... Link +/-45 plies together
%heq=[XDES(2)-XDES(3)];
heq=[XDES(1)-XDES(2)];

% Add constraint to history
%his_g_buckl=[his_g_buckl;g_bukl];
hisG=[hisG; max(max(g))];

% Plot results for this iteration
plotDesIter(his, hisG, wingDATA,lamDATA,solDATA,ldsDATA,desDATA,1);

% Plot thickness surface at iteration intervals of 10 function evaluations
if length(his)>=iPlotThkSurf*length(XDES)
    iPlotThkSurf=iPlotThkSurf+1;
    disp('PlotThkSurf = '); PlotThkSurf(iPlotThkSurf)
    if PlotThkSurf(iPlotThkSurf)==0

        plotTitle='Thickness Distribution';
        iterTitle=' FunEval=';  iterTitle=strcat(iterTitle,int2str(length(his)));

        Wpanel=panelWeightVthk(solDATA,lamDATA,wingDATA.panXY )
        weightTitle=' Weight=';
        weightTitle=strcat(weightTitle,num2str(Wpanel,6));
        
        plotTitle=strcat(plotTitle,',', iterTitle, ',', weightTitle);
        [psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(solDATA.Mg,solDATA.intMethod,wingDATA.panXY,0);
        getLayerThk(lamDATA.tDeg,psiI,etaI,lamDATA,1,plotTitle);
        %PlotThkSurf(iPlotThkSurf)=1;
    end        
end

%function [f,hx,xC]=vlam_Freq_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
function [f]=vlam_Freq_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Freq_Obj.m : plate frequency objective function

global his
global hisG

% original
%lamDATA.thk=XDES;
%lamDATA.thk_coeff=[1 0 0 0];

% lamDATA.thk=XDES(1:2);
% lamDATA.thk_coeff=[XDES(3) 0 0 0];
% lamDATA.thk_coeff=[XDES(1) 0 XDES(2) 0]; % t=linear with x
% lamDATA.thk_coeff=[XDES(1) XDES(3) XDES(2) XDES(4)]; % t=bi-linear with x & y

if (lamDATA.tDeg==0 && desDATA.dvType(1)==0)
    lamDATA.thk=XDES; % constant thickness
elseif (desDATA.dvType(1)==1)
    lamDATA.theta0(1)=XDES(1);
    lamDATA.theta0(2)=XDES(1)+45;
    lamDATA.theta0(3)=XDES(1)-45;
    lamDATA.theta0(4)=XDES(1)+90;
    lamDATA.theta1(1)=XDES(2);
    lamDATA.theta1(2)=XDES(2)+45;
    lamDATA.theta1(3)=XDES(2)-45;
    lamDATA.theta1(4)=XDES(2)+90;
elseif (lamDATA.tDeg>0 && desDATA.dvType(2)==2)
    lamDATA.thk_coeff=XDES; % variable thickness
end


% Eigenvalue solution
solDATA.isoltype=1;
[C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);

% Set objective function value: maximize 1st frequency
WN=sqrt(EgN);
f=-WN(1);

% Add current weight and thicknesses to design history his
%Wt = panelWeight(lamDATA.thk,lamDATA.rho_ply,wingDATA);
Wt = panelWeightVthk(solDATA,lamDATA,wingDATA.panXY );
%his=[his;XDES WN(1) Wt];
% his=[his;XDES(1:2) WN(1) Wt];
%his=[his;XDES WN(1)];
hisG=[hisG;XDES WN(1)];


function [f,hx,xC]=Buckl_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% Buckl_Obj.m : plate buckling eignevalue objective function

global his

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness

    % Fiber angles only    
    %setVarLamTheta( lamDATA, XDES(1), XDES(2) );    
    
    nply=length(lamDATA.theta0);
    for i=1:nply/2
        if(mod(i,2) == 1)
            lamDATA.theta0(i)=XDES(1);
            lamDATA.theta1(i)=XDES(2);
        else
            lamDATA.theta0(i)=-XDES(1);
            lamDATA.theta1(i)=-XDES(2);
        end
    end
    j=1;
    for i=nply:-1:nply/2+1
        lamDATA.theta0(i)=lamDATA.theta0(j);
        lamDATA.theta1(i)=lamDATA.theta1(j);
        j=j+1;
    end
else
    %lamDATA.thk_coeff=XDES; % variable thickness
end

% Eigenvalue solution
[C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);

% Set objective function value
f=-EgN(1);

% Add constraint to design history his
Wt = panelWeightVthk(solDATA,lamDATA,wingDATA.panXY );
[g,heq]=fibK_Con(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);

%his=[his;XDES EgN(1) Wt];
his=[his;XDES EgN(1) g];


function [f,hx,xC]=Strain_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% Strain_Obj.m : strain objective function

global his

if(lamDATA.tDeg==0)
    %lamDATA.thk=XDES; % constant thickness

    % Fiber angles only    
    %setVarLamTheta( lamDATA, XDES(1), XDES(2) );    
    
    nply=length(lamDATA.theta0);
    for i=1:nply/2
        if(mod(i,2) == 1)
            lamDATA.theta0(i)=XDES(1);
            lamDATA.theta1(i)=XDES(2);
        else
            lamDATA.theta0(i)=-XDES(1);
            lamDATA.theta1(i)=-XDES(2);
        end
    end
    j=1;
    for i=nply:-1:nply/2+1
        lamDATA.theta0(i)=lamDATA.theta0(j);
        lamDATA.theta1(i)=lamDATA.theta1(j);
        j=j+1;
    end
else
    %lamDATA.thk_coeff=XDES; % variable thickness
end


% Ritz solution
[C,ritzWt,EgN,ritzDisp,ritzEps,sMesh,tMesh]=vLam(wingDATA,lamDATA,solDATA,ldsDATA);

% Set objective function value: minimize max strain
ritzMax = [max(max(ritzEps(:,1))); max(max(ritzEps(:,2))); max(max(ritzEps(:,3))); max(max(ritzEps(:,4))); max(max(ritzEps(:,5)))];
ritzMin = [min(min(ritzEps(:,1))); min(min(ritzEps(:,2))); min(min(ritzEps(:,3))); min(min(ritzEps(:,4))); min(min(ritzEps(:,5)))];

% Get max strain
ritzMaxAbsAllComp = max(max(abs(ritzMax)));
ritzMinAbsAllComp = max(max(abs(ritzMin)));
ritzMaxAbs = max([ritzMaxAbsAllComp ritzMinAbsAllComp]);

f=ritzMaxAbs;

% Add constraint to design history his
Wt = panelWeightVthk(solDATA,lamDATA,wingDATA.panXY );
[g,heq]=fibK_Con(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);

%his=[his;XDES EgN(1) Wt];
his=[his;XDES f g];

function [g,heq]=strain_Con(Eps,wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA)
%% Strain_Obj.m : strain objective function

% Allowables
allow_epsx = .0018;
allow_epsy = .0018;
allow_epsxy= .0027;
allow_epsxz= .0015;
allow_epsyz= .0015;

% Find Min/Max strains
epsDim=size(Eps);
%     epsMax = [max(max(Eps(:,1))); max(max(Eps(:,2))); max(max(Eps(:,3))); max(max(Eps(:,4))); max(max(Eps(:,5)))];
%     epsMin = [min(min(Eps(:,1))); min(min(Eps(:,2))); min(min(Eps(:,3))); min(min(Eps(:,4))); min(min(Eps(:,5)))];
% limit to in-plane strains
epsMax = [max(max(Eps(:,1))); max(max(Eps(:,2))); max(max(Eps(:,3)))];
epsMin = [min(min(Eps(:,1))); min(min(Eps(:,2))); min(min(Eps(:,3)))];

% Get max strain
epsMaxAbsAllComp = max(max(abs(epsMax)));
epsMinAbsAllComp = max(max(abs(epsMin)));
epsMaxAbs = max([epsMaxAbsAllComp epsMinAbsAllComp]);

% Calc constraint
g=epsMaxAbs/0.001800 -1;

%%% Normalize constraints: Eps_act/Eps_all <= 1, or Eps_act/Eps_all-1 <=0

% Fiber Strains: Eps-x
g_epsx = abs(Eps(:,1)/allow_epsx) -1;

% Fiber Strains: Eps-y
g_epsy = abs(Eps(:,2)/allow_epsy) -1;

% Shear Strains: Eps-xy
g_epsxy = abs(Eps(:,3)/allow_epsxy) -1;

% Transverse Strain
g_epsxz = abs(Eps(:,4)/allow_epsxz) -1;

% Transverse Strain
g_epsyz = abs(Eps(:,5)/allow_epsyz) -1;

% Combine constraints
%g=[g_epsx; g_epsy; g_epsxy; g_epsxz; g_epsyz];
g=[g_epsx; g_epsy; g_epsxy];

% gmax=[max(max(g_epsx)); max(max(g_epsy)); max(max(g_epsxy))];
% gmin=[min(min(g_epsx)); min(min(g_epsy)); min(min(g_epsxy))];

% Use full size of strain vectors
%g=abs(Eps(:,2))/allow_epsy-1;


%disp('size of g_strain = '); size(g)

heq=0;  % Eq constraints



function [g,heq]=vlam_Volume_Con(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% vlam_Volume_Con.m - plate volume constraint function g

% Compute volume constraint value
% lamDATA.thk=XDES;
% lamDATA.thk=XDES(1:2);

if(lamDATA.tDeg==0 && desDATA.dvType(1)==0)
    lamDATA.thk=XDES; % constant thickness
elseif (lamDATA.tDeg>0 && desDATA.dvType(1)==0)
    lamDATA.thk_coeff=XDES; % variable thickness
end
%Wt = panelWeight(lamDATA.thk,lamDATA.rho_ply,panXY);
Wt = panelWeightVthk(solDATA,lamDATA,wingDATA.panXY );
WtBound = 2.0; % 1.0 lbs max weight
g=Wt/WtBound-1; %Ineq constraints

if(lamDATA.tDeg > 0)
    % Thickness constraints
    %[g_thk,heq_thk]=vlamThkCon(XDES,panXY,lamDATA,solDATA);
    
    % Combine panel weight constraint and thickness constraints into vector
    %g = [g; g_thk];
end

heq=0;  % Eq constraints



function [g,heq]=fibK_Con(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% fibK_Con.m - fiber curvature constraint function

panXY=wingDATA.panXY;

% Allowable curvature
Kallow=1/desDATA.Rallow;

% Start by checking only 1st ply
T0=XDES(1);
T1=XDES(2);

a=abs(panXY(2,1)-panXY(1,1));
x=linspace(0,a/2);

fibK = getFibCurvature(a,T0,T1,x);
fibK_min=min(fibK);

% Find curvature less than input allowable
I=find(fibK<=Kallow);

% Normalize constraint
%g=fibK/Kallow-1
%g=fibK_min/Kallow-1;
g=Kallow/fibK_min-1;

heq=0;  % Eq constraints


function [g,heq]=vlamThkCon(XDES,wingDATA,lamDATA,solDATA,meshDATA,intDATA)
%% vlamThkCon.m - plate thickness constrain function g

% Enforces postive thickness over plate surface

% Set lamDATA thk coeffs to optimized values
%lamDATA.thk_coeff=XDES;
lamDATA.thk_coeff=XDES(1,3:end);

if(lamDATA.tDeg > 0)
    % mesh data for thk evaluation pts
%     [ meshDATA, intDATA ] = meshGen( wingDATA, solDATA );
    [psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts,surfI,xyzI] = getIntData(intDATA);

    % compute surface representation of thickness
%     thk_coeff_mat = repmat(lamDATA.thk_coeff,20,1);
%     z1 = thk_legendre(lamDATA.tDeg,psiIMesh,etaIMesh)*thk_coeff_mat';
%     disp('size of lamDATA.thk_coeff='); size(lamDATA.thk_coeff)
%     disp('size of thk_coeff_mat='); size(thk_coeff_mat)
%     disp('size of z1='); size(z1)
%     disp('size of psiIMesh='); size(psiIMesh)
%     disp('size of etaIMesh='); size(etaIMesh)
%    surface(psiIMesh,etaIMesh,z1,'FaceColor',[0.5 1.0 0.5], 'EdgeColor', 'none');

%     dim = length(psiIMesh)^2;
%     psiIMesh2=reshape(psiIMesh,1,dim); etaIMesh2=reshape(etaIMesh,1,dim);
%     z1 = thk_legendre(lamDATA.tDeg,psiIMesh2,etaIMesh2)*lamDATA.thk_coeff';
%     z1 = reshape(z1,size(psiIMesh));
%     %surface(psiIMesh,etaIMesh,z1,'FaceColor',[0.5 1.0 0.5], 'EdgeColor', 'none');
%     
%     % define a plane at z=0
%     z2 = zeros(size(psiIMesh));
%     %surface(psiIMesh,etaIMesh,z2,'FaceColor',[1.0 0.5 0.0], 'EdgeColor', 'none');

    % Use the uniform mesh instead of integ pts
    [psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(20); % Increase point sampling to improve thickness constraint continuity
    dim = length(psiBCMesh)^2;
    psiBCMesh2=reshape(psiBCMesh,1,dim); etaBCMesh2=reshape(etaBCMesh,1,dim);
    z1 = thk_legendre(lamDATA.tDeg,psiBCMesh2,etaBCMesh2)*lamDATA.thk_coeff';
    z1 = reshape(z1,size(psiBCMesh));

    % define a plane at z=0
    z2 = zeros(size(psiBCMesh));
    
    % intersection of computed thk surface with min-thk surface
    zdiff = z1 - z2;
    zdiff_norm=zdiff/0.125 - 1.0; % normalize thickness (0.125 is arbitratry)
    %test_neg=all(zdiff < 0)
    
    % return values
    %g=reshape(zdiff_norm,length(zdiff)^2,1);
    tminGage=0.025;
    g=reshape(-zdiff/tminGage,length(zdiff)^2,1); 
    heq=0;
    
    % get min-max thicknesses
    tmax = max(max(zdiff));
    tmin = min(min(zdiff));

    % Draw thickness surface and contour of intersection
    if(solDATA.plot_flag > 0)
        figure
        surface(psiIMesh,etaIMesh,z1,'FaceColor',[0.5 1.0 0.5], 'EdgeColor', 'none');
        surface(psiIMesh,etaIMesh,z2,'FaceColor',[1.0 0.5 0.0], 'EdgeColor', 'none');
        
        %  Contours of surface intersection
        C = contours(psiIMesh, etaIMesh, zdiff, [0 0]);
        %disp('Size of contour intersections = '); C
        %disp('Norm of C = ');norm(C)
        
        % Extract the x- and y-locations from the contour matrix C.
        xL = C(1, 2:end);
        yL = C(2, 2:end);
        
        % Interpolate on the first surface to find z-locs for the intersection
        zL = interp2(psiIMesh, etaIMesh, z1, xL, yL);
        
        % Draw the line of intersection between surfaces.
        line(xL, yL, zL, 'Color', 'k', 'LineWidth', 3);     
        
        title('\bfSkin Thickness Distribution')
        ylabel('\bfThickness')
        
        view( 127.5, 30 );

    end
end

function [A]=getInequalLayerCon(wingDATA,lamDATA,solDATA,meshDATA,intDATA)
%% getInequalLayerCon.m - plate thickness constrain function g

% Inequality constraint Matrix A - enforces postive thickness over plate surface
if(lamDATA.tDeg > 0)
    % mesh data for thk evaluation pts
    %[ meshDATA, intDATA ] = meshGen( wingDATA, solDATA );
    [psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts,surfI,xyzI] = getIntData(intDATA);

%     dim = length(psiIMesh)^2;
%     psiIMesh2=reshape(psiIMesh,1,dim); etaIMesh2=reshape(etaIMesh,1,dim);
%     A = thk_legendre(lamDATA.tDeg,psiIMesh2,etaIMesh2)';

    %[psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(10);
    [psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(20); % Increase point sampling to improve thickness constraint continuity
    dim = length(psiBCMesh)^2;
    psiBCMesh2=reshape(psiBCMesh,1,dim); etaBCMesh2=reshape(etaBCMesh,1,dim);
    A = thk_legendre(lamDATA.tDeg,psiBCMesh2,etaBCMesh2)';
end



function [XOPT,his,hisF]= Weight_Eign_Des_Opt(wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% Weight_Des_Opt.m : Min weight design optimization w/ strain constraint
disp('Min weight design & max freq optimization w/ strain constraint');

global his
his=[];

global hisF
hisF=[];

global hisG
hisG=[];

% Handles to objective and constraint functions
OBJ=@vib_weight_obj; % must return a vector
CON=@genCon;

wfactor=[1;1];
goal=[1.0;100];

% Thickness inequality constraints
A=[];b=[];
if( lamDATA.tDeg > 0 )
    [A]=getInequalLayerCon(wingDATA,lamDATA,solDATA)
     tmax(size(A))= 0.150;
     tmin(size(A))= 0.001;

%   Min thk constraint
    A=-A';
    b=-tmin;
    
%   Max & Min thk constraint
%     A=[A';-A'];
%     b=[tmax;-tmin];
end

% Recover design data from structure
[XL,XU,X0]=getDesData(desDATA);
disp('XL, XU, X0 ='); XL, XU, X0

% Optimize
OP=[];
OP=optimset('disp','iter')  % Display optimization/iteration details
%XOPT=fgoalattain(OBJ,X0,goal,wfactor,A,b,[],[],XL,XU,CON,OP,panXY,lamDATA,solDATA,ldsDATA,desDATA)
XOPT=fgoalattain(OBJ,X0,goal,wfactor,A,b,[],[],XL,XU,CON,OP,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)

% Place design and constraint history into one matrix
%his= [his hisG];

%function [weight_obj, vib_freq_obj] = vib_weight_obj(XDES,panXY,lamDATA,solDATA,ldsDATA,desDATA)
function [out] = vib_weight_obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA)
%% vibration and weight objective function
weight_obj=Weight_Des_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
vib_freq_obj=vlam_Freq_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);

out=[weight_obj; vib_freq_obj];


function plotDesIter=plotDesIter(his,hisG, wingDATA,lamDATA,solDATA,ldsDATA,desDATA,plot_hist)
%% Plot design solutions

global h_histFig;

% his = history matrix (design variables, obj fn, panel weight)
%       his(1:2) = design variables x(1), x(2)
%       his(3)   = obj fn
%       his(4)   = panel weight

%format long; may need for complex eigenvalue
format short;
% HIS=his
% display most current history values
%HIS=[his(end,:) hisG(end,:)]

numX=length(desDATA.X0);

% Determine if plot should be updated
dim=size(his);
nrows=dim(1);

%if(plot_hist > 0 && ~rem(nrows,numX)*nrows/numX > 0)
if(plot_hist > 0)
    
    % xvals = integers representing iteration
    xvals=linspace(1,nrows,nrows);

    ip=1;
    nvar=0;
    nP=length(desDATA.dvType)+2;
    
    % Check if figure is already created
    %figure
    if h_histFig == 0
        h_histFig = figure;
    else
       set(0, 'CurrentFigure', h_histFig);
       first_plot=1;
    end

    % Design variable history
    for i=1:length(desDATA.dvType)
        
        % Thickness variables
        if desDATA.dvType(i)==0
            subplot(nP,1,ip)
            nthk=length(lamDATA.thk);
            
 %           if first_plot == 0 
                plot(his(:,nvar+1:nvar+nthk),'linewidth',3)
                title('\bfDesign history')
                ylabel('\bfLayer Thks')
%            else
%                 set(h_histFig,'XData',xvals);
%                 set(h_histFig,'YData',his(:,nvar:nvar+nthk-1));
%                 refreshdata
%                 drawnow
%             end
            ip=ip+1;
            nvar=nvar+nthk;
        end
        
        % Fiber angle variables
        if desDATA.dvType(i)==1
            subplot(nP,1,ip)
            ntheta=2; % hardwired
            
%            if first_plot == 0
                plot(his(:,nvar+1:nvar+ntheta),'linewidth',3)
                ylabel('\bfTheta DVs')
%            else
%                 set(h_histFig,'XData',xvals);
%                 set(h_histFig,'YData',his(:,nvar+1:nvar+ntheta));
%                 refreshdata
%                 drawnow
%             end
            ip=ip+1;
            nvar=nvar+ntheta;
        end
        
        % Polynomial coefficients
        if desDATA.dvType(i)==2
            subplot(nP,1,ip)
            npoly=length(lamDATA.thk_coeff);
            
            %if first_plot == 0
                plot(his(:,nvar+1:nvar+npoly),'linewidth',3)
                ylabel('\bfPolyCoeff DVs')
            %else
%                 set(h_histFig,'XData',xvals);
%                 set(h_histFig,'YData',his(:,nvar:nvar+npoly-1));        
%                 refreshdata
%                 drawnow
            %end
            ip=ip+1;
            nvar=nvar+npoly;
        end
    end
    
    % Objective funtion history
    subplot(nP,1,ip)
    plot(his(:,nvar+1),'linewidth',3)
    %hold on;
    ylabel('\bfObjective')
    ip=ip+1;
    
    % Constraint function history
    subplot(nP,1,ip)
    plot(hisG(:,:),'linewidth',3)
    ylabel('\bfConstraint')
    xlabel('\bfFunction evaluation number')

    print('design_history', '-dpng');
       
end


function plotDes=plotDes(XOPT,his,wingDATA,lamDATA,solDATA,ldsDATA,desDATA,plot_hist,istart)
%% Plot design solutions

disp('Plotting design history');

% his = history matrix (design variables, obj fn, panel weight)
%       his(1:2) = design variables x(1), x(2)
%       his(3)   = obj fn
%       his(4)   = panel weight

format long;
HIS=his

numX=length(XOPT);

% if(plot_hist > 0)
%     figure
%     % Design variable history
%     subplot(3,1,1)
%     %plot(his(:,1:2),'linewidth',3)
%     %plot(his(:,1:4),'linewidth',3)
%     plot(his(:,1:numX),'linewidth',3)
%     title('\bfDesign history')
%     ylabel('\bfDesign variables')
%     
%     % Objective funtion history
%     subplot(3,1,2)
%     %plot(his(:,3),'linewidth',3)
%     %plot(his(:,5),'linewidth',3)
%     plot(his(:,numX+1),'linewidth',3)
%     hold on;
%     if(strcmp(desDATA.drespType,'WEIGHT_EIGN'))
%         plot(his(:,numX+2),'linewidth',3)
%         hold off;
%         
%         % Pareto Front
%         subplot(3,1,3)
%         plot( his(:,numX+1), his(:,numX+2),'-o','linewidth',2 )
%         
%         if(lamDATA.tDeg > 0)
%             solDATA.plot_flag=1;
%             [g,heq]=vlamThkCon(XOPT,wingDATA,lamDATA,solDATA);
%         end
% 
%         return;
%     end
%     ylabel('\bfObjective')
%     
%     % Constraint function history
%     subplot(3,1,3)
%     %plot(his(:,4),'linewidth',3)
%     %plot(his(:,6),'linewidth',3)
%     size(his)
%     plot(his(:,numX+2),'linewidth',3)
%     %plot(his(:,numX+1),'linewidth',3)  % Freq Problem
%     ylabel('\bfConstraint')
%     xlabel('\bfFunction evaluation number')
% 
%     print('design_history', '-dpng');
%end


if(plot_hist > 0)

    ip=1;
    nvar=0;
    nP=length(desDATA.dvType)+2;
    
    figure
    
    % Design variable history
    for i=1:length(desDATA.dvType)
        
        % Thickness variables
        if desDATA.dvType(i)==0
            subplot(nP,1,ip)
            nthk=length(lamDATA.thk);
            plot(his(:,nvar+1:nvar+nthk),'linewidth',3)
            title('\bfDesign history')
            ylabel('\bfLayer Thks')
            ip=ip+1;
            nvar=nvar+nthk;
        end
        
        % Fiber angle variables
        if desDATA.dvType(i)==1
            subplot(nP,1,ip)
            ntheta=2; % hardwired
            plot(his(:,nvar+1:nvar+2),'linewidth',3)
            ylabel('\bfTheta DVs')
            ip=ip+1;
            nvar=nvar+ntheta;
        end
        
        % Polynomial coefficients
        if desDATA.dvType(i)==2
            subplot(nP,1,ip)
            npoly=length(lamDATA.thk_coeff);
            plot(his(:,nvar+1:nvar+npoly),'linewidth',3)
            ylabel('\bfPolyCoeff DVs')
            ip=ip+1;
            nvar=nvar+npoly;
        end
    end
    
    % Objective funtion history
    subplot(nP,1,ip)
    plot(his(:,nvar+1),'linewidth',3)
    hold on;
    ylabel('\bfObjective')
    ip=ip+1;
    
    % Constraint function history
    subplot(nP,1,ip)
    plot(his(:,nvar+2),'linewidth',3)
    ylabel('\bfConstraint')
    xlabel('\bfFunction evaluation number')

    print('design_history', '-dpng');        
       
end

% Compute data for objective function contours
global his
his=[];

global hisG
hisG=[];

% Range of design variables for Obj Fn contours
x1=linspace(desDATA.XL(1),desDATA.XU(1)+0.1,10);
x2=linspace(desDATA.XL(2),desDATA.XU(2)+0.1,10);

% Early attempt to use polynomial coeffs
% x1=linspace(-1.0,1.0,10);
% x2=linspace(-1.0,1.0,10);
% x3=linspace(-1.0,1.0,10);
% x4=linspace(-1.0,1.0,10);


%% Contour plot of objective function and constraints
if(lamDATA.tDeg==0 && numX==2)
    if(istart ==1)  % only do this once if mult. start pts
        [xx,yy]=meshgrid(x1,x2);
        for i=1:length(x1)
            for j=1:length(x2)
                XDES=[x1(i) x2(j)]; %original
                
                %XDES=[x1(i)  x2(j) thk_coeff(1)];
                %XDES=[x1(i)  x2(j)]; % linear with x
                %XDES=[x1(i)  x2(j)  x3(j)  x4(j)];  % bi-linear with x,y
                %XDES=[x1(i)  0.00  0.00  x2(j)  0.00  0.00   x3(j)  0.00  0.00];  % bi-linear with x,y
                
                % compute obj fn value
                switch desDATA.drespType
                    case {'WEIGHT'}
                        f=Weight_Des_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
                        g=genCon(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
                    case {'EIGN'}
                        f=vlam_Freq_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
                    case {'BUCKL'}
                        f=-Buckl_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
                    case {'STRAIN'}
                        f=Strain_Obj(XDES,wingDATA,lamDATA,solDATA,ldsDATA,desDATA);
                    otherwise
                        f=0;
                        return;
                end
            end
        end
        his = [his hisG];
        CHIS=his;  % Note: 'his' is updated automatically from Obj Fn calls above
        
        % Objective Fn Contours: OBJC
        figure
        OBJC=reshape(his(:,numX+1),size(xx));  %his(:,numX+1) = obj fn
        objmin=min(min(OBJC));objmax=max(max(OBJC));
        [cc,hh]=contour(xx,yy,OBJC);  % transposed OBJC
        set(hh,'linewidth',3)
        clabel(cc)
        hold on
        
        %     % ========================================================
        %     % Draw feasible design region based on min steering radius
        %     % ========================================================
        %     a=abs(panXY(2,1)-panXY(1,1));
        %     x=linspace(0,a/2);
        %     xx1=linspace(desDATA.XL(1),desDATA.XU(1),25);
        %     xx2=linspace(desDATA.XL(2),desDATA.XU(2),25);
        %     [xxc,yyc]=meshgrid(xx1,xx2);
        %
        %     for i=1:length(xx1)
        %         for j=1:length(xx2)
        %             fibK = getFibCurvature(a,xx1(i),xx2(j),x);
        %             fibK_min(i,j) = min(min(abs(fibK)));
        %             fibK_max(i,j) = max(max(abs(fibK)));
        %         end
        %     end
        %     % Find curvatures less than input allowable
        %     I=find(fibK_min<=1/desDATA.Rallow);
        %
        %     % Plot allowable fiber curvature points
        %     X=xxc(:);Y=yyc(:);
        %     plot(X(I),Y(I),'y.','MarkerSize',50);
        %     % ========================================================
        
        % Constraint Contours: CONC
        %CONC=reshape(his(:,numX+2),size(xx))  %his(:,numX+2) = constraint
        CONC=reshape(his(:,numX+1),size(xx))  %his(:,numX+2) = constraint
        gmax=max(max(CONC))
        gmin=min(min(CONC))
        if(gmax ~= gmin)
            if(gmax > 10)
                gmax=1;
            end
            CL=linspace(gmin,gmax,10);
            CL = [CL 0]; % ensure that pt for feasible constraint boundary included in level
            CL = sort(CL);
            %CL=linspace(gmin,gmax,5);
            [cc1,hh1]=contour(xx,yy,CONC,CL);clabel(cc1)    % all contours
            
            %[cc1,hh1]=contour(xx,yy,CONC,[0 0]);clabel(cc1)  % constraint boundary = [0 0]
            set(hh1,'linewidth',2,'color','k')
            
        else
        end
        
        %     % Redraw obj contours to overlay on feasible region
        %     [cc,hh]=contour(xx,yy,OBJC);  % transposed OBJC
        %     set(hh,'linewidth',3)
        %     clabel(cc)
        
    else
        hold on
    end
    
    % Overlay optimized result values
    plot(XOPT(1),XOPT(2),'rp','linewidth',5)
    plot(HIS(:,1),HIS(:,2),'b--','linewidth',3)
    
    % plot(XOPT(1),XOPT(4),'rp','linewidth',5)
    % plot(HIS(:,1),HIS(:,4),'b--','linewidth',3)
    
    % title({'\bfObjective Function and Constraint Contours'; ...
    %        '\bfTwo-Layer Cantilevered Plate'; ...
    %        '\bfAluminum/CarbonEpoxy'});
    [panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(wingDATA.panXY, solDATA, lamDATA);
    
    title({'\bfObjective Function and Constraint Contours'; ...
        panTitle; lamTitle});
    xlabel('\bfDesign Variable X_1');
    ylabel('\bfDesign Variable X_2');
    
    hold off
    print('design_space_coutours', '-dpng');
    
else
    %Draw panel thickness function
    plotTitle='';
    [psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(solDATA.pDeg,solDATA.intMethod,wingDATA.panXY,solDATA.plot_flag);
    %lamDATA.thk_coeff=XOPT;
    lamDATA.thk_coeff=XOPT(1,3:end);
    layer_thk = getLayerThk(lamDATA.tDeg,psiI,etaI,lamDATA,1,plotTitle);
end

function stop = outfun(x,optimValues,state)
%% output function
     stop = false;
 
     switch state
         case 'init'
             hold on
         case 'iter'
         % Concatenate current point and objective function
         % value with history. x must be a row vector.
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x; x];
         % Concatenate current search direction with 
         % searchdir.
           searchdir = [searchdir;... 
                        optimValues.searchdirection'];
           %plot(x(1),x(2),'o');
         % Label points with iteration number and add title.
         % Add .15 to x(1) to separate label from plotted 'o'
           %text(x(1)+.15,x(2),... 
           %     num2str(optimValues.iteration));
           %title('Sequence of Points Computed by fmincon');
         case 'done'
             hold off
         otherwise
     end
