function [solDATA,lamDATA,panXY,ldsDATA] = readProblem(problem,solDATA,lamDATA,panXY,ldsDATA)

switch problem
    case 1  % Vib frequency convergence: QI laminate, SSSS
        paramType = 1; 
        solDATA.soltype=1;
        bcType='ssss';
        lamDATA.ilam=2;
end