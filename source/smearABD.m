function [ ABD_smear ] = smearABD( ABD, thk )
%smearABD calculates smeared laminate ABD properties
%   Detailed explanation goes here

[Ex Ey nuxy Gxy] = lam_engr_constants(abd, thk);
Exy = Ex
ABD_smear = ABD;

t = sum(thk);

Dx = t^3*Ex/(1-.;
Dy = t^3*Ey/12.;
Dxy= t^3*

end

