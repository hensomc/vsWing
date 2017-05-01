function [ desDATA ] = setDesData(XL,XU,X0,Xstart,Rallow,drespType,dvType)
%setDesData Sets design obtmization data structure

%  desDATA = returned data structure

desDATA.XL = XL;
desDATA.XU = XU;
desDATA.X0 = X0;
desDATA.Xstart = Xstart;

desDATA.Rallow = Rallow;

desDATA.drespType = drespType;
desDATA.dvType = dvType;