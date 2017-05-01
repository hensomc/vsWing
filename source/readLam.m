function [imat, theta0, theta1, thk, mat_names, phi_rot] = readLam(ilam)
%function [imat, theta0, theta1, thk, mat_names]

% purpose: read laminate stack data into memory

    %disp('readLam: ilam='); ilam
    
    % default
    phi_rot = 0.;
    
    switch ilam
        case 0
            % Case 1: 1-ply isotropic plate
            mat_names={'aluminum'};
            imat=[1];
            theta0=[0];
            %thk=[.1];
            thk=[0.2];
            %thk=[0.025];
            theta1=[0];
        
        case 1
            % Case 1: 2-ply isotropic plate
            mat_names={'aluminum'};
            imat=[1 1];
            theta0=[0 0];
            thk=[.05 .05];
            theta1=[0 0];
            
        case 2
            % Case 2: 8-ply Quasi layup
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[0 45 -45 90 90 -45 45 0];
            %theta0=[0 0 0 0 0 0 0 0];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053 ];
            %theta1=[0 0 0 0 0 0 0 0];
            theta1=[0 45 -45 90 90 -45 45 0];
            
        case 3
            % Case 3: 2-ply curvilinear layup
            mat_names={'carbon_ep'};
            T0=45;
            T1=0;
            imat=[1 1];
            theta0=[T0 -T0];
            thk=[.0053 .0053];
            theta1=[T1 -T1];
            
        case 4
            % Case 4: 8-ply curvilinear layup
            mat_names={'carbon_ep'};
            %mat_names={'hi_ortho_carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[ 45 -45  45 -45 -45  45 -45  45];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053];
            theta1=[ 0 0 0 0 0 0 0 0];
            
        case 5
            % Case 5: 4-ply unsymmetic layup, see Houmat 2013 pp211-224
            mat_names={'carbon_ep'};
            imat=[1 1 1 1];
            theta0=[+45 -45 +45 -45];
            thk=[.01 .01 .01 .01];
            theta1=[+45 -45 +45 -45];
            
        case 6
            % Case 6: 2-ply bi-metallic laminate
            mat_names={'aluminum', 'steel'};
            imat=[1 2];
            theta0=[0 0];
            thk=[0.05 0.05];
            theta1=[0 0];
            
        case 7
            % Case 7: 2-ply aluminum-carbon_ep laminate
            mat_names={'aluminum', 'carbon_ep'};
            imat=[1 2];
            theta0=[0 0];
            thk=[0.05 0.05];
            theta1=[0 0];
            
        case 8
            % Case 8: 8-ply Uni layup
            %mat_names={'carbon_ep'};
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[0 0 0 0 0 0 0 0];
            %theta0=[0 0 0 0 0 0 0 0];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053 ];
            %theta1=[0 0 0 0 0 0 0 0];
            theta1=[0 0 0 0 0 0 0 0];
            
        case 9
            % Case 9: 8-ply +/- layup
            %mat_names={'carbon_ep'};
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[45 -45 45 -45 -45 45 -45 45];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053 ];
            theta1=[45 -45 45 -45 -45 45 -45 45];
            
        case 10
            % Case 10: 32-ply Quasi layup
            mat_names={'carbon_ep'};
            imat(1:32)=1;
            theta0=[45 -45 0 90 45 -45 0 90 45 -45 0 90 45 -45 0 90 ...
                    90 0 -45 45 90 0 -45 45 90 0 -45 45 90 0 -45 45];
            thk(1:32)=.0053;
            theta1=[45 -45 0 90 45 -45 0 90 45 -45 0 90 45 -45 0 90 ...
                    90 0 -45 45 90 0 -45 45 90 0 -45 45 90 0 -45 45];
                
        case 11
            % Case 11: 2-ply +/-45 layup
            mat_names={'carbon_ep'};
            imat=[1 1];
            theta0=[45 -45];
            thk=[.0053 .0053];
            theta1=[45 -45];
            
        case 12
            % Case 12: 2-ply +/-30 layup
            mat_names={'hi_ortho_carbon_ep'};
            imat=[1 1];
            theta0=[30 -30];
            thk=[.053 .053];
            theta1=[30 -30];
            
        case 13
            % Case 13: 8 ply curvilinear, %[<-45|45>,<45|-45>,<-45|45>,<45|-45>,  <45|-45>,<-45|45>,<45|-45>,<-45|45>]
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[-45 +45 -45 +45  +45 -45 +45 -45];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053];
            theta1=[+45 -45 +45 -45  -45 +45 -45 +45];

        case 14
            % Case 14: Gurdal 2008-Case 1 laminate: 8 ply curvilinear, %[<0|45>,<0|-45>,<0|45>,<0|-45>,  <0|-45>,<0|45>,<0|-45>,<0|45>]
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[0 0 0 0 0 0 0 0];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053];
            theta1=[+45 -45 +45 -45  -45 +45 -45 +45];       
            phi_rot = 30;
            
        case 15
            % Case 15: Gurdal 2008-Case 1 laminate: 8 ply curvilinear, %[<90|45>,<90|-45>,<90|45>,<90|-45>,  <90|-45>,<90|45>,<90|-45>,<90|45>]
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[90 90 90 90 90 90 90 90];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053];
            theta1=[+45 -45 +45 -45  -45 +45 -45 +45];
            
        case 16
            % Case 16: 32-ply Angle Ply layup
            mat_names={'carbon_ep'};
            imat(1:32)=1;
            theta0=[45 -45 45 -45 45 -45 45 -45 45 -45 45 -45 45 -45 45 -45 ...
                    -45 45 -45 45 -45 45 -45 45 -45 45 -45 45 -45 45 -45 45 ];
            thk(1:32)=.0053;
            theta1=[45 -45 0 90 45 -45 0 90 45 -45 0 90 45 -45 0 90 ...
                90 0 -45 45 90 0 -45 45 90 0 -45 45 90 0 -45 45];
            
        case 17
            % Case 17: 8 ply curvilinear, %[<-45|45>,<45|-45>,<-45|45>,<45|-45>,  <45|-45>,<-45|45>,<45|-45>,<-45|45>]
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[-45 +45 -45 +45  +45 -45 +45 -45];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053];
            theta1=[0 0 0 0 0 0 0 0];
            phi_rot = 70;
            
        case 20
            % Case 20: sandwich panel with Quasi facesheets
            mat_names={'carbon_ep';'syncore'};
            imat=[1 1 1 1 2 2 1 1 1 1];
            theta0=[45 -45 0 90 0 0 90 0 -45 45];
            thk=[.0053 .0053 .0053 .0053 0.50 0.50 .0053 .0053 .0053 .0053];
            theta1=[45 -45 0 90 0 0 90 0 -45 45];
            phi_rot = 0;
            
        case 21
            % Case 21: 8 ply curvilinear, wing skewed
            mat_names={'carbon_ep'};
            imat=[1 1 1 1 1 1 1 1];
            theta0=[-45 +45 -45 +45  +45 -45 +45 -45];
            thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053];
            theta1=[20 20 20 20 20 20 20 20];
            %phi_rot = 45;
            
        case 22
            % Case 22: 4-ply curvilinear smeared layup, variable 0-degree
            mat_names={'carbon_ep'};
            T0=0;
            T1=-30;
%             T0=30;
%             T1= 0;
            imat=[1 1 1 1];
%             theta0=[T0 45 -45 90 90 -45 45 T0];
%             theta1=[T1 60 -60 120 120 -60 60 T1];
%             thk=[.0053 .0053 .0053 .0053 .0053 .0053 .0053 .0053 ];
%             thk=[.053 .053 .053 .053 .053 .053 .053 .053 ];
            theta0=[T0 T0+45 T0-45 T0+90];
            %theta1=[T1 60 -60 120];
            theta1=[T1 T1+45 T1-45 T1+90];
            %thk=[.0053 .0053 .0053 .0053];
            %thk=[.053 .053 .053 .053];
            %thk=[.005 .005 .005 .005];
            thk=[.1 .1 .1 .1];
%             dx=.0001;
%             thk=[.1+sqrt(-1)*dx .1 .1 .1];
            phi_rot = 0;

        case 23
            % Case 23: 4-Ply Quasi for Smeared Laminate
            mat_names={'carbon_ep'};
            imat=[1 1 1 1];
            plyThk=0.0053;
            nPlyLay=10;
            tLay=nPlyLay*plyThk;
            theta0=[90 90 90 90];
            theta1=[90 90 90 90];
            thk=[tLay tLay tLay tLay];
            phi_rot = 0;
        
        case 106
            % Case 1: 1-ply isotropic plate
            mat_names={'aluminum'};
            imat=[1];
            theta0=[0];
            thk=[1.8];
            theta1=theta0;

            
    end
end