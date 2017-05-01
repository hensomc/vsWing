function [degenGeom,panXY] = importDegenGeom( dgfile, plot_flag )
%importDgenGeom Reads OpenVSP degengeom file
%   Detailed explanation goes here

  [degenGeom,panXY] = vspPlotDegenStick(dgfile, plot_flag);
  % hardwire panXY to panXY(1)
  panXY=cell2mat(panXY(1));
  
  if plot_flag>0
      vspPlotDegenPlate(dgfile);
      vspPlotDegenSurf(dgfile);
  end
end

