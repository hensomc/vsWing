function [C, EgN]=ritzSolModes( solDATA, Kq, Mq, T, plot_flag )

if(plot_flag > 0)
    disp('Solving for eigenvalues...');
end

bcMethod=solDATA.bcMethod;

% Solve
[PHI,EG]=eig(Kq,Mq); % ----------------------------------------- Eq(26)
[Eg,Is]=sort(diag(EG));
Phi=PHI(:,Is);
EgN=sqrt(Eg);

% Compute Ritz polynomial coeffs for each mode
C=cell(length(EgN));
for im=1:length(EgN)
    if(bcMethod == 1)
        C{im}=Phi(:,im);
    elseif(bcMethod == 2)
        C{im}=T*Phi(:,im);
    end
end