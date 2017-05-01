function [ wtvals, eigvals, dispvec, epsvec, jobStatus ] = readNasF06( solDATA, f06_filename )
%readNasF06 Reads Nastran F06 file into memory using system call to Perl
%   Detailed explanation goes here


%% Recover solution parameters
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);


%% Default return arguments
eigvals = 0;
wtvals = 0;
numGrid=(Mg+1)^2;
numElem=Mg^2;
dispvec = zeros(numGrid,6);
epsvec = zeros(numGrid,6);
e_strnvec = zeros(numElem,16);
g_strnvec = zeros(numGrid,16);
jobStatus = 0;


%% Scan f06 file
if(plot_flag > 0)
    disp('scanning f06 file for results');
end

fid = fopen(f06_filename);

ifound_mass_wt=0;
ifound_ev=0;
ifound_disp=0;
ifound_strn=0;
ifound_anal_summary=0;

tline = fgets(fid);

% Results counters
iwt=0;
iev=0;
idisp=0;
ielem=0;

% Read to end of file
while ischar(tline)
    %disp(tline)
    % Check for FATALs
    if( strfind(tline, 'FATAL') )
        disp('FATAL error found in f06 file:');tline
        jobStatus = -1;
        return;
    end

    % Mass Weight Check
    if( ifound_mass_wt == 0 )
        if(strfind(tline, 'MASS AXIS SYSTEM (S)     MASS              X-C.G.        Y-C.G.        Z-C.G.'));
            ifound_mass_wt=1;
            if(plot_flag > 0)
                disp('found mass weight check header in f06');
            end
        end
    end    
    
    % Eigenvalue header
    if( ifound_ev == 0 )
        if (strfind(tline, 'NO.       ORDER' ));
            ifound_ev=1;
            if(plot_flag > 0)
                disp('found eigenvalue marker in f06');
            end
        end
    end
    
    % Displacement header
    if( ifound_disp == 0 )
        if(strfind(tline, 'POINT ID.   TYPE          T1             T2             T3             R1             R2             R3'));
            ifound_disp=1;
            if(plot_flag > 0)
                disp('found displacement header in f06');
            end
        end
    end
    
    % Strain header
    if( ifound_strn == 0 )
        if(strfind(tline, 'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )'));
            ifound_strn=1;
            if(plot_flag > 0)
                disp('found strain header in f06');
            end
        end
    end    
    
    % read next line
    tline = fgets(fid);
    %disp(tline)

    % Break and return if numMass found
    if(iwt == 3) % 3 Rows of output in F06
        ifound_mass_wt = 0;
        %break;
    end
    
    % Break and return if numModes found
    if(iev == solDATA.numMode)
        ifound_ev = 0;
        %break;
    end
    
    % Break and return if numGrid found
    if(idisp == numGrid)
        ifound_disp = 0;
        %break;
    end

    % Break and return if numElem found
    if(ielem == numElem)
        ifound_strn = 0;
    end

    % Read mass weight summary
    if( ifound_mass_wt == 1 )
        iwt=iwt+1;
        vals = [sscanf(tline, '%s%f%f%f%f')];
        %disp({'vals[1:5] = ',vals});
        wtvals(iwt,1:4)= vals(2:5);
    end
    
    % Read eigenvalues
    if( ifound_ev == 1 )
        iev=iev+1;
        vals = [sscanf(tline, '%i%i%f%f%f%f%f')];
        %disp({'vals[3] = ',vals(3)});
        eigvals(iev,1)= sqrt(vals(3));
    end

    % Read displacements
    if( ifound_disp == 1 )
        C = strsplit(tline);
        if(length(C) == 10)
            idisp=idisp+1;
            vals = sscanf(tline, '%i%s%f%f%f%f%f%f');
            dispvec(idisp,1:6) = vals(3:8);
        end
    end
    
    % Analysis summary header
    if( ifound_anal_summary == 0 )
        if(strfind(tline, 'A N A L Y S I S  S U M M A R Y  T A B L E'));
            ifound_anal_summary=1;
            if(plot_flag > 0)
                disp('found analysis summary in f06');
            end
        end
    end    
    
    % Read strains
    if( ifound_strn == 1 && ifound_anal_summary == 0)
        C = strsplit(tline);
        if(length(C) == 10 || length(C) == 11 || length(C) == 12) 
            %print('tline='),tline

            % Skip miscellaneous lines
            if(strcmp(C(1),'1'))
                
            elseif(length(C) == 12 && strcmp(C(1),'0') )
                %vals = sscanf(tline, '%i%i%s%f%f%f%f%f%f%f%f');
                ielem=ielem+1;
                %e_strnvec(ielem,1:8) = vals(4:11)
                for i=1:8
                    e_strnvec(ielem,i)=str2num(C{1,i+3});
                end
                e_strnvec(ielem,:);
                %pause
                next_rec = 'elem';
            
            
            elseif(length(C) == 10 && strcmp(next_rec,'elem'))
                %vals = sscanf(tline, '%s%f%f%f%f%f%f%f%f');                
                %e_strnvec(ielem,9:16) = vals(2:9)
                for i=1:8
                    e_strnvec(ielem,i+9)=str2num(C{1,i+1});
                end
                %e_strnvec(ielem,:);
                %pause
                next_rec = 'grid';
            
            
            elseif(length(C) == 11 && strcmp(next_rec,'grid'))
                %vals = sscanf(tline, '%s%i%f%f%f%f%f%f%f%f');
                gid=str2num(C{1,2});
                for i=1:8
                    vals2(1,i)=str2num(C{1,i+2});
                end
                if(max(max(g_strnvec(gid,1:8)))==0.0)
                    a=vals2(1,1:8);
                else
                    a=mean([g_strnvec(gid,1:8);vals2(1,1:8)]);
                end
                g_strnvec(gid,1:8) = a(1:8);
                
                %pause
                next_rec = 'grid';
               
            
            elseif(length(C) == 10 && strcmp(next_rec,'grid'))
                %vals = sscanf(tline, '%s%f%f%f%f%f%f%f%f');  
                for i=1:8
                    vals3(1,i)=str2num(C{1,i+1});
                end
                if(max(max(g_strnvec(gid,9:16)))==0.0)
                    a=vals3(1,1:8);
                else
                    a=mean([g_strnvec(gid,9:16); vals3(1,1:8)]);
                end
                
                g_strnvec(gid,9:16) = a(1:8);
                
                %pause
                next_rec = 'elem';
            end      
            
        end
    end
end

% assign Nodal strain vector: epsvec
%epsvec(:,1:3)=g_strnvec(:,10:12);  % Assign only top fiber strain for now
epsvec(:,1:6)=g_strnvec(:,10:15);  % Assign only top fiber strain for now, include prin strains & theta
fclose(fid);

if(plot_flag>0)
    disp('Finished reading f06');
end

end

