function [ meshDATA, intDATA ] = meshGen( wingDATA, solDATA )
%meshGen Generates mesh of points for Ritz plate analyses
%   Detailed explanation goes here

% Recover solution parameters
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

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

%% Uniform BC mesh in psi-eta domain
[psiBCMesh,etaBCMesh,isoGridID]=bcMeshIso(Mg);
[xBCMesh,yBCMesh] = meshTransform(panXY, psiBCMesh, etaBCMesh);

%% Integration point mesh
[psiIMesh,etaIMesh,psiI,etaI,gl_wts]=getIntPts(pDeg,intMethod,panXY,plot_flag);

%% Integration point mesh projected to upper/lower wing surfaces
if isempty(degenGeom)
    xyzInt=zeros(length(psiI),3);
    xyzI(1).s=struct('xyz',xyzInt);
    surfI=[];
else
    for i=1:nSurf
        surfI(i)=getDegenWingSurf_fvc(degenGeom,i,plot_flag);
        xyzInt=getWingZCoords(surfI(i),psiI,etaI,panXY,i-1,plot_flag);  % integ pts on wing surface
        xyzI(i).s=struct('xyz',xyzInt);
    end
end

%% Return meshDATA
meshDATA.psiBCMesh = psiBCMesh;
meshDATA.etaBCMesh = etaBCMesh;
meshDATA.isoGridID = isoGridID;

meshDATA.xBCMesh = xBCMesh;
meshDATA.yBCMesh = yBCMesh;

intDATA.psiIMesh = psiIMesh;
intDATA.etaIMesh = etaIMesh;
intDATA.psiI = psiI;
intDATA.etaI = etaI;
% intDATA.xIMesh = xIMesh;
% intDATA.yIMesh = yIMesh;
intDATA.gl_wts = gl_wts;

intDATA.surfI = surfI;
intDATA.xyzI = xyzI;

end