function setVarLamTheta( lamDATA, T0, T1)
% Set layup angles <T0|T1> based on inputs
% Laminate symmetry enforced

nply=length(lamDATA.theta0);
for i=1:nply/2
    if(mod(i,2) == 1)
        lamDATA.theta0(i)=T0;
        lamDATA.theta1(i)=T1;
    else
        lamDATA.theta0(i)=-T0;
        lamDATA.theta1(i)=-T1;
    end
end
j=1;
for i=nply:-1:nply/2+1
    lamDATA.theta0(i)=lamDATA.theta0(j);
    lamDATA.theta1(i)=lamDATA.theta1(j);
    j=j+1;
end
% 
% lamDATA.theta0
% lamDATA.theta1
% pause

