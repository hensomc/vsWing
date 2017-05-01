function [psiIMesh,etaIMesh,psiI,etaI,gl_wts,surfI,xyzI] = getIntData(intDATA)

psiIMesh = intDATA.psiIMesh;
etaIMesh = intDATA.etaIMesh;
psiI = intDATA.psiI; 
etaI = intDATA.etaI; 
% xIMesh = intDATA.xIMesh; 
% yIMesh = intDATA.yIMesh; 
gl_wts = intDATA.gl_wts; 

surfI =intDATA.surfI;
xyzI = intDATA.xyzI;

end