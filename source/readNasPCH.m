%function readNasPCH(pchFile);
%clear all

%fid=fopen('your_file.pch','r'); 
%fid=fopen(pchFile,'r'); 
fid=fopen('plate_fem.pch','r'); 
block_size = 100;

a=1;
while ~feof(fid) 
          tline = textscan(fid,'%s',block_size,'delimiter','\n');
          e = 1;
          while e<=length(tline{1,1})
               l = length(tline{1,1}{e,1}); 
               for k = 1:(ceil(l/8))
                b = k*8;
                if (l)>=b
                  PCH{a,1}{e,k} = tline{1,1}{e,1}((b-7):b);
                else 
                  PCH{a,1}{e,k} = tline{1,1}{e,1}((b-7):(l));
                end
               end
               e = e+1;
           end
    a = a+1; 
end
fclose(fid);