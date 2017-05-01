function[edgeBC,edgeGID]=setEdgeBC(Mg, bcType)
%% Sets panel edge boundary conditions

% Edge Grid IDs
edge1GID=[1:Mg+1];
edge3GID=[Mg*(Mg+1)+1 : Mg*(Mg+1)+(Mg+1)];
for(i=1:Mg+1);
   edge2GID(i) = (Mg+1)*(i-1) +1;
   edge4GID(i) = (Mg+1)*i;
end
edgeGID=[edge1GID; edge2GID; edge3GID; edge4GID];

% Decode BCs
for(i=1:length(bcType))
    if(bcType(i) == 's')  % simple
        edgeBC(i,:) = [1 0 0];
    elseif(bcType(i) == 'c')  % clamped
        edgeBC(i,:) = [1 1 1];
    elseif(bcType(i) == 'f')  % free
        edgeBC(i,:) = [0 0 0];
    end
end

% Corner BCs
% c1=1; c2=Mg*(Mg+1)+1; c3=(Mg+1)^2; c4=Mg+1; 
% cornerGID=[c1 c2 c3 c4];
% cornerBC=[];  % make a 4x3 [u v w; u v w; u v w; u v w]