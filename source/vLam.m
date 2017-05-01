function [C,ritzWt,EgN,ritzDisp,ritzEps]=vLam(wingDATA,lamDATA,solDATA,ldsDATA,meshDATA,intDATA)
%vLam - main driver for variable stiffness laminate

% start timer
tic

% default return values
EgN=0;Kq=0;Mq=0;PHI=0;EG=0;
ritzDisp=0;ritzEps=0;C=0;
Fc=0;

% Recover panel data
panXY = wingDATA.panXY;
degenGeom = wingDATA.degenGeom;

% Surfaces
nsurf=0;
if isempty(degenGeom)
    nSurf=1;
else
    nSurf=2;
end

% Echo laminate data
if(solDATA.plot_flag > 0)
    lamDATA
end

%% Recover solution parameters
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);
parforArg=Inf;
%parforArg=0;

%% Recover mesh data
[psiBCMesh,etaBCMesh,isoGridID,xBCMesh,yBCMesh] = getMeshData(meshDATA);
[psiIMesh,etaIMesh,psiI,etaI,gl_wts,surfI,xyzI] = getIntData(intDATA);

%% Set Edge Boundary Condition Grid IDs
[edgeBC,edgeGID]=setEdgeBC(Mg, bcType);

%% PRE solution plots
preSolPlots(plot_flag, wingDATA, lamDATA, solDATA, psiBCMesh, etaBCMesh, xBCMesh, yBCMesh, psiI, etaI);

%%  Generate basis functions
[B]=getRitzBasisVectors(ipoltype,pDeg,psiI,etaI,plot_flag); % ---------- Eq(13)

%% Plate Stiffness and Mass Matrices
ndof=5*(pDeg+1)^2;
Kplate=zeros(ndof,ndof);Mplate=zeros(ndof,ndof);
Wplate=0.0;
if isempty(degenGeom)
    zI=zeros(length(psiI),3);
    [Kplate,Mplate,Wplate]=KMplate(pDeg,intMethod,smearKM,panXY,lamDATA,psiI,etaI,zI,gl_wts,B,plot_flag);
    surfI=[];
else
    %for i=1:nSurf  % later - consider generalize on 1 surface sweep of lifting surface
    %parfor i=1:nSurf, parforArg; % later - consider generalize on 1 surface sweep of lifting surface
    for i=1:nSurf % later - consider generalize on 1 surface sweep of lifting surface
        xyzInt=xyzI(i).s.xyz;
        zIskin(i,:)=xyzInt(:,3);     
        [K,M,Wp]=KMplate(pDeg,intMethod,smearKM,panXY,lamDATA,psiI,etaI,xyzInt(:,3)',gl_wts,B,plot_flag);
        
        Kplate=Kplate+K;
        Mplate=Mplate+M;
        Wplate=Wplate+Wp;
    end
end


%KM_matrix_qualities(Kplate,Mplate,2)

%% Spar cap stiffness & mass matrices
Kscap=zeros(ndof,ndof); Mscap=zeros(ndof,ndof);
Wtcap=0.0;
zIweb=[];

% spars
sz=size(wingDATA.spar.XY);
nSpar=sz(2);
%if nSpar>0 && wingDATA.sparCapXsect>0
if nSpar>0
    sparCapLamDATA = setLamData(wingDATA.sparCapLamID, 0, []);

    bmW=2; % beamWidth - move to sparDATA section variable
    %eqW=[bmW/(panXY(2,1)-panXY(1,1)) bmW/(panXY(3,1)-panXY(4,1)) 0];
    eqW=[bmW bmW 0];
    
    
    %for j=1:nSurf  % upper/lower surfaces
    %parfor j=1:nSurf, parforArg;  % upper/lower surfaces
    parfor j=1:nSurf;  % upper/lower surfaces
        for i=1:nSpar
            beamXY=cell2mat(wingDATA.spar.XY(i));
            [psiIMeshBeam,etaIMeshBeam,psiIbeam,etaIbeam,gl_wts_beam]=getIntPtsBeam(Mg,intMethod,beamXY,panXY,plot_flag);
            if nSurf==1
                xyzIbeam=zeros(length(etaI),3);
            else
                xyzIbeam=getWingZCoords(surfI(j),psiIbeam,etaIbeam,panXY,j-1,plot_flag);
            end
            
            if wingDATA.sparCapXsect>0
                %plot3(xyzIbeam(:,1),xyzIbeam(:,2),xyzIbeam(:,3),'o')
                [Bbeam]=getRitzBasisVectors(ipoltype,pDeg,psiIbeam,etaIbeam,plot_flag);
                [Kb,Mb,Wb]=KMbeam(pDeg,intMethod,smearKM,panXY,sparCapLamDATA,psiIbeam,etaIbeam,xyzIbeam(:,3)',gl_wts_beam,Bbeam,eqW,plot_flag);
                Kscap = Kscap + Kb;
                Mscap = Mscap + Mb;
                Wtcap = Wtcap + Wb;
            end
            zIweb(i,j,:)=xyzIbeam(:,3)';
        end
    end
    if nSurf == 1
        zIweb(:,2,:)=zIweb(:,1,:);
    end
end

%% Spar web stiffness & mass matrices
Ksweb=zeros(ndof,ndof); Msweb=zeros(ndof,ndof);
Wtweb=0.0;

if nSpar>0 && wingDATA.sparWebXsect>0
%     if isempty(degenGeom)
%     else
        sparWebLamDATA = setLamData(wingDATA.sparWebLamID, 0, []);
        webThk=sum(sparWebLamDATA.thk); % web thickness - move to sparDATA section variable
        %webThk=sum(sparWebLamDATA.thk); % web thickness - move to sparDATA section variable
        %eqT=[webThk webThk 0];  % effective web thickness
        %eqT=[webThk/(panXY(2,1)-panXY(1,1)) webThk/(panXY(3,1)-panXY(4,1)) 0];  % effective web thickness
        eqT=[webThk webThk 0];  % effective web thickness

        %for i=1:nSpar
        %parfor i=1:nSpar,parforArg;
        parfor i=1:nSpar;
            beamXY=cell2mat(wingDATA.spar.XY(i));
            [psiIMeshBeam,etaIMeshBeam,psiIbeam,etaIbeam,gl_wts_beam]=getIntPtsBeam(Mg,intMethod,beamXY,panXY,plot_flag);
            [Bbeam]=getRitzBasisVectors(ipoltype,pDeg,psiIbeam,etaIbeam,plot_flag);
            % Offsets for spar webs and no upper/lower surfaces
%             if isempty(zIweb)
%                 zIweb(1:nSpar,1,1:length(psiI))=1.0;
%                 zIweb(1:nSpar,2,1:length(psiI))=-1.0;
%             end
            [Kb,Mb,Wb]=KMweb(pDeg,intMethod,smearKM,panXY,sparWebLamDATA,psiIbeam,etaIbeam,zIweb(i,1,:),zIweb(i,2,:),gl_wts_beam,Bbeam,eqT,plot_flag);
            Ksweb = Ksweb + Kb;
            Msweb = Msweb + Mb;
            Wtweb = Wtweb + Wb;
        end
%     end
end
%disp('Max stiffness terms for [plate, cap, web]: '); max(max(Kplate)), max(max(Kscap)), max(max(Ksweb))

%% Rib Webs

sz=size(wingDATA.rib.XY);
nRib=sz(2);

if nRib>0
    disp('nRibs = '); nRib
end

%% Bay Materials
Kbay=zeros(ndof,ndof); Mbay=zeros(ndof,ndof); Wtbay=0;
if length(wingDATA.bayMatName) > 0
    [Kbay,Mbay,Wtbay]=KMbay(pDeg,intMethod,panXY,lamDATA,psiI,etaI,zIskin(1,:)',zIskin(2,:)',gl_wts,B,{wingDATA.bayMatName},plot_flag);
end

%% Sum all stiffnesses & masses
Kc=Kplate+Kscap+Ksweb+Kbay;
Mc=Mplate+Mscap+Msweb+Mbay;


%% Collect Weight Data
ritzWt.Wplate=Wplate;
ritzWt.Wtbay=Wtbay;
ritzWt.Wtcap=Wtcap;
ritzWt.Wtweb=Wtweb;
ritzWt.WtTot=Wplate+Wtbay+Wtcap+Wtweb;


%% Potential energy due to external loads (geometric stiffness matrix)
if(isoltype == 2 || isoltype == 3)
    [G]=K_ExtLds(panXY,isoltype,ipoltype,pDeg,intMethod,psiI,etaI,gl_wts,B,ldsDATA,plot_flag);
    if(isoltype == 2)
        Mc=-G;
    elseif(isoltype == 3)
        Fc=G;
    end
end

%% Apply BCs
[Kq,Mq,Fq,T]=applyBC( panXY,solDATA,edgeBC,edgeGID,isoGridID,bcType,ldsDATA,Kc,Mc,Fc,plot_flag );

%% Solve eigenvalue problem
if(isoltype == 1 || isoltype == 2)
    [C, EgN]=ritzSolModes( solDATA, Kq, Mq, T, plot_flag );

% Solve linear statics problem
elseif( isoltype == 3 )  % Solve [Kq]{q}={Fq}  --------------------- Eq(35)
    [C, ritzDisp, ritzEps]=ritzSolLinear( solDATA, panXY, lamDATA, psiBCMesh(:), etaBCMesh(:), surfI, Kq, Fq, T, plot_flag );
end

%% Stop timer
toc