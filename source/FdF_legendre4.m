function [F,Fx,Fy,Fxy,Fxx,Fyy] = FdF_legendre4(nDeg,x,y)

F   = F_legendre_4(x,y);
Fx  = Fx_legendre_4(x,y);
Fy  = Fy_legendre_4(x,y);
Fxy = Fxy_legendre_4(x,y);
Fxx = Fxx_legendre_4(x,y);
Fyy = Fyy_legendre_4(x,y);


function F = F_legendre_4(x,y)
%F_LEGENDRE_4
%    F = F_LEGENDRE_4(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:22:06

t2 = y.^2;
t3 = t2.*(3.0./2.0);
t4 = t3-1.0./2.0;
t5 = t2.*5.0;
t6 = t5-3.0;
t7 = t2.^2;
t8 = t7.*(3.5e1./8.0);
t13 = t2.*(1.5e1./4.0);
t9 = t8-t13+3.0./8.0;
t10 = x.^2;
t11 = t10.*(3.0./2.0);
t12 = t11-1.0./2.0;
t14 = t10.*5.0;
t15 = t14-3.0;
t16 = t10.^2;
t17 = t16.*(3.5e1./8.0);
t19 = t10.*(1.5e1./4.0);
t18 = t17-t19+3.0./8.0;
F = [1.0,y,t4,t6.*y.*(1.0./2.0),t9,x,x.*y,t4.*x,t6.*x.*y.*(1.0./2.0),t9.*x,t12,t12.*y,t4.*t12,t6.*t12.*y.*(1.0./2.0),t9.*t12,t15.*x.*(1.0./2.0),t15.*x.*y.*(1.0./2.0),t4.*t15.*x.*(1.0./2.0),t6.*t15.*x.*y.*(1.0./4.0),t9.*t15.*x.*(1.0./2.0),t18,t18.*y,t4.*t18,t6.*t18.*y.*(1.0./2.0),t9.*t18];


function Fx = Fx_legendre_4(x,y)
%FX_LEGENDRE_4
%    FX = FX_LEGENDRE_4(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:22:18

t2 = y.^2;
t3 = t2.*(3.0./2.0);
t4 = t3-1.0./2.0;
t5 = t2.*5.0;
t6 = t5-3.0;
t7 = t2.^2;
t8 = t7.*(3.5e1./8.0);
t13 = t2.*(1.5e1./4.0);
t9 = t8-t13+3.0./8.0;
t10 = x.^2;
t11 = t10.*(1.5e1./2.0);
t12 = t11-3.0./2.0;
t14 = t10.*x.*(3.5e1./2.0);
Fx = [0.0,0.0,0.0,0.0,0.0,1.0,y,t4,t6.*y.*(1.0./2.0),t9,x.*3.0,x.*y.*3.0,t4.*x.*3.0,t6.*x.*y.*(3.0./2.0),t9.*x.*3.0,t12,t12.*y,t4.*t12,t6.*t12.*y.*(1.0./2.0),t9.*t12,t14-x.*(1.5e1./2.0),y.*(t14-x.*(1.5e1./2.0)),t4.*(t14-x.*(1.5e1./2.0)),t6.*y.*(t14-x.*(1.5e1./2.0)).*(1.0./2.0),t9.*(t14-x.*(1.5e1./2.0))];


function Fy = Fy_legendre_4(x,y)
%FY_LEGENDRE_4
%    FY = FY_LEGENDRE_4(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:22:47

t2 = y.^2;
t3 = t2.*(1.5e1./2.0);
t4 = t3-3.0./2.0;
t5 = t2.*y.*(3.5e1./2.0);
t6 = x.^2;
t7 = t6.*(3.0./2.0);
t8 = t7-1.0./2.0;
t9 = t6.*5.0;
t10 = t9-3.0;
t11 = t6.^2;
t12 = t11.*(3.5e1./8.0);
t14 = t6.*(1.5e1./4.0);
t13 = t12-t14+3.0./8.0;
Fy = [0.0,1.0,y.*3.0,t4,t5-y.*(1.5e1./2.0),0.0,x,x.*y.*3.0,t4.*x,x.*(t5-y.*(1.5e1./2.0)),0.0,t8,t8.*y.*3.0,t4.*t8,t8.*(t5-y.*(1.5e1./2.0)),0.0,t10.*x.*(1.0./2.0),t10.*x.*y.*(3.0./2.0),t4.*t10.*x.*(1.0./2.0),t10.*x.*(t5-y.*(1.5e1./2.0)).*(1.0./2.0),0.0,t13,t13.*y.*3.0,t4.*t13,t13.*(t5-y.*(1.5e1./2.0))];


function Fxy = Fxy_legendre_4(x,y)
%FXY_LEGENDRE_4
%    FXY = FXY_LEGENDRE_4(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:22:54

t2 = y.^2;
t3 = t2.*(1.5e1./2.0);
t4 = t3-3.0./2.0;
t5 = t2.*y.*(3.5e1./2.0);
t6 = x.^2;
t7 = t6.*(1.5e1./2.0);
t8 = t7-3.0./2.0;
t9 = t6.*x.*(3.5e1./2.0);
Fxy = [0.0,0.0,0.0,0.0,0.0,0.0,1.0,y.*3.0,t4,t5-y.*(1.5e1./2.0),0.0,x.*3.0,x.*y.*9.0,t4.*x.*3.0,x.*(t5-y.*(1.5e1./2.0)).*3.0,0.0,t8,t8.*y.*3.0,t4.*t8,t8.*(t5-y.*(1.5e1./2.0)),0.0,t9-x.*(1.5e1./2.0),y.*(t9-x.*(1.5e1./2.0)).*3.0,t4.*(t9-x.*(1.5e1./2.0)),(t9-x.*(1.5e1./2.0)).*(t5-y.*(1.5e1./2.0))];


function Fxx = Fxx_legendre_4(x,y)
%FXX_LEGENDRE_4
%    FXX = FXX_LEGENDRE_4(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:23:05

t2 = y.^2;
t3 = t2.*5.0;
t4 = t3-3.0;
t5 = t2.^2;
t6 = x.^2;
t7 = t6.*(1.05e2./2.0);
t8 = t7-1.5e1./2.0;
t9 = t2.*(3.0./2.0);
t10 = t9-1.0./2.0;
t11 = t5.*(3.5e1./8.0);
t12 = t2.*(-1.5e1./4.0)+t11+3.0./8.0;
Fxx = [0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,y.*3.0,t2.*(9.0./2.0)-3.0./2.0,t4.*y.*(3.0./2.0),t2.*(-4.5e1./4.0)+t5.*(1.05e2./8.0)+9.0./8.0,x.*1.5e1,x.*y.*1.5e1,t10.*x.*1.5e1,t4.*x.*y.*(1.5e1./2.0),t12.*x.*1.5e1,t8,t8.*y,t8.*t10,t4.*t8.*y.*(1.0./2.0),t8.*t12];

function Fyy = Fyy_legendre_4(x,y)
%FYY_LEGENDRE_4
%    FYY = FYY_LEGENDRE_4(X,Y)

%    This function was generated by the Symbolic Math Toolbox version 6.0.
%    25-Apr-2015 21:23:20

t2 = y.^2;
t3 = t2.*(1.05e2./2.0);
t4 = t3-1.5e1./2.0;
t5 = x.^2;
t6 = t5.*(3.0./2.0);
t7 = t6-1.0./2.0;
t8 = t5.*5.0;
t9 = t8-3.0;
t10 = t5.^2;
t11 = t10.*(3.5e1./8.0);
t12 = t5.*(-1.5e1./4.0)+t11+3.0./8.0;
Fyy = [0.0,0.0,3.0,y.*1.5e1,t4,0.0,0.0,x.*3.0,x.*y.*1.5e1,t4.*x,0.0,0.0,t5.*(9.0./2.0)-3.0./2.0,t7.*y.*1.5e1,t4.*t7,0.0,0.0,t9.*x.*(3.0./2.0),t9.*x.*y.*(1.5e1./2.0),t4.*t9.*x.*(1.0./2.0),0.0,0.0,t5.*(-4.5e1./4.0)+t10.*(1.05e2./8.0)+9.0./8.0,t12.*y.*1.5e1,t4.*t12];



