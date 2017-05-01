function [ solDATA ] = setSolData( ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType )

%setSolData Sets the solution parameters data structure
%  ipoltype      % 1=power,  2=legendre, 3=power-all numerical
%  pDegI;        % polynomial degre
%  plot_flag     % 0 = none, 1=modes, 2=animate modes, 3=details
%  isoltype      % 1=vib,    2=buckl,    3=linear static
%  iDesOpt       % 0 = analysis only, 1= frequency design optimization
%  intMethod     % 1=simpson rule, 2=gauss quadrature
%  bcMethod      % 1=spring, 2=explicit null space, 3=spring general
%  bcType='abcd' % 4 edge BC's of panel, 'c'=clamped, 's'=simple, 'f'=free
%  RitzSE        % 0=default=off, 1=on
%  numMode       % 15 = # of modes to recover
%  ssCalc        % 0 = skip, 1 = compute stress/strain recovery
%  paramType     % 0 = none, 1 = Poly Order, 2= Integ Pts Mg, 3 = Plate Skew, 4 = Layup Angle <T0,T1>

%  solDATA = returned data structure

solDATA.ipoltype=ipoltype;
solDATA.pDeg=pDeg;
solDATA.Mg=Mg;
solDATA.smearKM=smearKM;
solDATA.plot_flag=plot_flag;
solDATA.isoltype=isoltype;
solDATA.iDesOpt=iDesOpt;
solDATA.femSol=femSol;
solDATA.symSol=symSol;
solDATA.intMethod=intMethod;
solDATA.bcMethod=bcMethod;
solDATA.bcType=bcType;
solDATA.RitzSE=RitzSE;
solDATA.bcSpring=bcSpring;
solDATA.iThkFun=iThkFun;
solDATA.numMode=numMode;
solDATA.ssCalc=ssCalc;
solDATA.paramType=paramType;