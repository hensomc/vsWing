%% Revise Mg to build FE mesh consistently for Spars & Ribs
function Mg = check_Mg(nSpar,Mg_i)

valid=0;
Mg=Mg_i;

if nSpar ==0
    return
end

% Mg/(nSpar+1) must be int: mod(Mg,nSpar)==0
nnx=Mg/(nSpar+1)+1;

% Case: Odd # spars
if(mod(nSpar,2)==1)
    if(mod(Mg,nSpar+1)==0 && nnx >= 2)
        return
    else
        while valid == 0
            Mg=Mg+1;
            if(mod(Mg,nSpar+1)==0)
                valid=1;
            elseif(Mg>=100)
                break;
            end
        end
    end
% Case: even # spars    
else
    if(mod(Mg,nSpar+1)==0 && nnx >= 2 && mod(nnx,2)==1)
        return
    else
        while valid == 0
            Mg=Mg+1;
            nnx=Mg/(nSpar+1)+1;
            if(mod(Mg,nSpar+1)==0 && mod(nnx,2)==1)
                valid=1;
            elseif(Mg>=100)
                break;
            end
        end
    end
end