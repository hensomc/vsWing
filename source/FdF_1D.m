function [F,dFx,dFxx]=FdF_1D(ipoltype,Nmax,xn)
%% Purpose:  Generate F and its derivative, numerically, vectorized operation

if(ipoltype == 1)
    [F,dFx,dFxx]=FdF_power_1D(Nmax,xn);
end

if(ipoltype == 2)
    [F,dFx,dFxx]=FdF_legendre_1D(Nmax,xn);
end

end