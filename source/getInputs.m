function [wingDATA,solDATA,lamDATA,ldsDATA,desDATA]=getInputs(jobID,degenGeomFileName)
%% Get inputs

% solver defaults
[solDATA]=setSolDefaults();
[ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc,paramType]=getSolData(solDATA);

% wing geometry defaults
degenGeom=[];
rsLayout = 0;
nRib = 0;
nSpar = 0;
sparCapXsect = 0;
sparWebXsect = 0;
sparWebLamID = 0;
sparCapLamID = 0;
bayMatName = '';

% default loads vectors
ldsDATA=setLdsDefaults();
[ Nx, Ny, Nxy, ptLds, ptLoc, p0 ] = getLdsData( ldsDATA );

% lam data default - constant thickness
tDeg=0; % constant thickness

% design variable defaults
[XL,XU,X0,Xstart,drespType,Rallow,dvType]=setDesDefaults();

% Job inputs
switch jobID
    case 0 % degenGeomFileName provided
        panID=0;
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
        ilam=22;
        bcType='fcff';
        disp('degenGeomFileName = '); degenGeomFileName
        
    case 1 % Param Study: Free Vib QI: Freq Convergence (pDeg: 4-12)
        iDesOpt=2;
        isoltype = 1;
        paramType=1;
        Mg = 20;
        femSol=0;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'cccc'; 'cfff';
        ilam=2;
        
    case 2 % Param Study: Free Vib VAT1: Freq Integ Order Convergence (Mg: 10-40)
        iDesOpt=2;
        isoltype=1;
        pDeg = 10;
        femSol=0;
        paramType=2;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'cccc'; 'cfff';
        ilam=13;
        
    case 3 % Analysis: Free Vib VAT1: Ritz to FEM comparison
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        pDeg = 12;
        Mg = 20;    

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'cccc'; 'cfff';
        ilam=17; % 13, 14, 15, 17, 2, 9
        plot_flag=1;
        
    case 4 % Analysis: Buckl: VAT1 Ritz to FEM comparison - Axial compression
        iDesOpt=0;
        isoltype=2;
        femSol=1;
        pDeg = 8;  % Note: use pDeg=8 for ssss and ilam=2 - otherwise complex e.v.
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cscs'; %'ssss'; 'cccc'; 'cfff'; 'cfcf';
        ilam=17; % 13, 14, 15, 2, 9
        %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
        numMode=5;
        plot_flag=1;
        Nx=  0; Ny=-100; Nxy=0; % Ny buckling case

    case 5 % Analysis: Static VAT1: Ritz to FEM comparison - Cantilevered plate with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'cfff'; %'ssss'; 'cccc';
        ilam=13; %13;
        p0=-0.025; % uniform pressure
        plot_flag=1;
        
    case 6 % Analysis: Static VAT1: Ritz to FEM comparison - In plane bending - tip shear
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'ssss'; 'cccc';
        ilam=13; %13, 15
        plot_flag = 1;
        Nx=0; Ny=0; Nxy=100; % tip shear case

    case 7 % Param Study: Free Vib: Skew Plate Geom (skew = 0-60)
        iDesOpt=2;
        isoltype=1;
        paramType=3;
        femSol=1;
        pDeg = 8;
        Mg = 20;
        numMode=5;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'ssss' 'cccc'; 'cfff';
        ilam=13;  % 13
       
    case 8 % Analysis: Free Vib: VAT2 Quad Plate Geom
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        pDeg = 12;
        Mg = 20;

        panID=10;
        panXY = [-4.0 -4.0;  % Quadrilateral, same as above, shifted
            6.0 -4.0;
            4.0  6.0;
            -4.0  4.0];

        bcType='cfff'; %'cccc'; 'cfff';
        ilam=15; % 13, 14, 15, 9, 2

    case 9 % Analysis: Buckling: VAT2 Quad Plate Geom
        iDesOpt=0;
        isoltype=2;
        femSol=1;
        pDeg = 8;
        Mg = 20;
        numMode = 5;

        panID=10;
        panXY = [-4.0 -4.0;  % Quadrilateral, same as above, shifted
            6.0 -4.0;
            4.0  6.0;
            -4.0  4.0];

        bcType='cccc'; %'ssss'; 'cccc'; 'cfff';
        ilam=9;  % 13, 14, 15, 2, 9
        %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
        Nx=   0; Ny=-10; Nxy=0; % Ny buckling case

        
    case 10 % Analysis: Static: VAT2 Quad Plate Geom
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 10;
        Mg = 20; %60

        panID=10;
        panXY = [-4.0 -4.0;  % Quadrilateral, same as above, shifted
            6.0 -4.0;
            4.0  6.0;
            -4.0  4.0];

        bcType='cccc'; %'ssss; 'cccc'; 'cfff';
        ilam=9; % 13, 14, 15, 2, 9
        p0=-0.025;
        plot_flag=0;
        
%     case 11 % Analysis: Static: VAT2 Quad Plate Geom
%         iDesOpt=0;
%         isoltype=3;
%         femSol=1;
%         pDeg = 8;
%         Mg = 20;
%         panID=4;
%         panXY = [ 0.0  0.0;  % Square
%             10.0  0.0;
%             10.0 10.0;
%             0.0 10.0];
%         bcType='ssss'; %'ssss'; 'cccc'; 'cfff';
%         ilam=9;  % 13, 14, 15, 2, 9
%         %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
%         Nx=   0; Ny=-100; Nxy=0; % Ny buckling case
%         plot_flag=1;

    case 11 % Param Study: Buckling VAT1: (T0:T1)
        iDesOpt=2;
        isoltype=2;
        pDeg = 6;
        femSol=0;
        paramType=4;
        numMode = 1;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'cccc'; 'cfff';
        ilam=9;
        Nx=-10; Ny=   0; Nxy=0; % Nx buckling case
        %Nx=   0; Ny=-10; Nxy=0; % Ny buckling case       
 
    case 12 % Optimization: Nx-Buckling: Des Vars T0 & T1
        iDesOpt=1;
        isoltype=2;
        femSol=0;
        pDeg = 3;
        Mg = 20;
        numMode = 5;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'ssss'; 'cccc'; 'cfff';
        ilam=9;  % 13, 14, 15, 2, 9
        %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
        Nx=  -10; Ny=0; Nxy=0; % Ny buckling case
        
    case 13 % Optimization: Ny-Buckling: Des Vars T0 & T1
        iDesOpt=1;
        isoltype=2;
        femSol=0;
        pDeg = 3;
        Mg = 20;
        numMode = 5;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cccc'; %'ssss'; 'cccc'; 'cfff';
        ilam=9;  % 13, 14, 15, 2, 9
        %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
        Nx=   0; Ny=-10; Nxy=0; % Ny buckling case

    case 14 % Buckling Test for Speed Up Work: Des Vars T0 & T1
        iDesOpt=0;
        isoltype=2;
        femSol=0;
        pDeg = 6;
        Mg = 20;
        numMode = 5;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'ssss'; 'cccc'; 'cfff';
        ilam=13;  % 13, 14, 15, 2, 9
        %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
        Nx=  -10; Ny=0; Nxy=0; % Ny buckling case
        plot_flag = 0;

        case 15 % Analysis: Static VAT1: Ritz to FEM comparison - Cantilevered plate with tip loads
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];
 
        bcType='cfff'; 
        ilam=1; %13;
        ptLds = [   0.   0.   -0.5;
                    0.   0.   -0.5 ];
        ptLoc = [  10.   0.;
                   10.  10.];
        plot_flag=1;
               
     case 16 % Analysis: Free Vib: Isotropic Quad Plate Geom - For DrWang
%         iDesOpt=0;
%         isoltype=1;
%         femSol=1;
%         pDeg = 8;
%         Mg = 20;
%         panID=10;
%         bcType='ssss'; %'cccc'; 'cfff';
%         ilam=0; %
%         plot_flag=1;
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        pDeg = 8;
        Mg = 20;

        panID=13;
        panXY = [ 0.0  0.0;  % Dr Wang test plate
            30.0  0.0;
            30.0 20.0;
            0.0  20.0];    

        bcType='ssss'; %'cccc'; 'cfff';
        ilam=0; %
        plot_flag=1;

    case 20 % Analysis: Static Sandwich Panel: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        numMode = 15;
        pDeg = 6;
        Mg = 10;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'cfff'; %'ssss'; 'cccc';
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=20; % sandwich panel
        
%         p0=-0.025; % uniform pressure
        
        % Tip torque
        ptLds = [   0.   0.   -0.5;
                    0.   0.    0.5 ];
        ptLoc = [   0.   10.;
                    10.  10.];
        
        plot_flag=1;
        
     case 21 % Analysis: Drawing fiber paths: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 20;

        panID=4; %square
        panID=8; %X-bias rect
        panID=14; %y-bias rect

        bcType=  'fcff'; %'cfff' ;'fcff'; 'ssss'; 'cccc';
        %ilam=0; % 
        ilam=13; % 
        ilam=22; % 
        p0=-0.005; % uniform pressure
        plot_flag=1; 

    case 96 % Cantilevered Unswept Wing y-dir span
        iDesOpt=1;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 10;
        
        panID=99;
        panXY = [ 0.0  0.0;  % Unswept rectangular wing y-dir
            50.0  0.0;
            50.0 110.0;
            0.0  110.0];
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alumn
        p0=-0.25; % uniform pressure
        plot_flag=1;

                % Design data
        drespType = 'WEIGHT';
        % ... Thickness, Linear in y
%         tDeg=1;
%         XL=[0.02  -1.0  0.0  0.0];
%         XU=[1.50   1.0  0.0  0.0];

        % ... Thickness, Quadratic in y    
        %F = [1.0, y, y^2, x, x.*y ,y^2*x ,x^2 ,x^2*y ,y^2*x^2];
%         tDeg=2;
%         XL=[ 0.005 -1.00  -1.00  0.00  0.00  0.00   0.00   0.00  0.00];
%         XU=[ 1.50   1.00   1.00  0.00  0.00  0.00   0.00   0.00  0.00];

        % ... Thickness, Cubic in y    
%         tDeg=3;
%         XL=[ 0.005 -1.00  -1.00  -1.00  0.00  0.00   0.00   0.00  0.00  0.00  0.00   0.00   0.00  0.00 0.00  0.00];
%         XU=[ 1.50   1.00   1.00   1.00  0.00  0.00   0.00   0.00  0.00  0.00  0.00   0.00   0.00  0.00 0.00  0.00];

        % ... Thickness, Quartic in y
%         tDeg=4;
%         XL=[ 0.005 -1.00  -1.00  -1.00 -1.00  -1.00  -1.00  -1.00 -1.00  -1.00  -1.00  -1.00   -1.00  -1.00 -1.00 -1.00 -1.00  -1.00 -1.00  -1.00   -1.00   -1.00  -1.00 -1.00  -1.00];
%         XU=[ 1.50   1.00   1.00   1.00  1.00  1.00   1.00   1.00  1.00  1.00  1.00   1.00   1.00  1.00 1.00  1.00  1.00  1.00  1.00  1.00   1.00   1.00  1.00 1.00  1.00];

        % .. Thickness, quadratic in x & y
        tDeg=2;
        XL=[ 0.005 -1.00  -1.00  -1.00  -1.00  -1.00   -1.00   -1.00  -1.00];
        XU=[ 0.75   1.00   1.00   1.00   1.00   1.00    1.00    1.00   1.00];

        % ... Thickness, Cubic in x & y    
%         tDeg=3;
%         XL=[ 0.005 -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00];
%         XU=[ 1.50   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00];
        
        X0=0.5*(XL+XU);
        thk_coeff=X0;
        
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs
        
     case 97 % Analysis: Cantilevered Unswept Wing y-dir span: SPARS WEBS at Midplane, Ritz to FEM comparison 
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 10;
        
        panID=99;
        panXY = [ 0.0  0.0;  % Unswept rectangular wing y-dir
            50.0  0.0;
            50.0 110.0;
            0.0  110.0];
        
        rsLayout=1234;
        nRib=0;
        nSpar=3;
        sparWebXsect=1;
        sparWebLamID=0;
        sparCapLamID=0;
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alumn
        %ilam=20; % sandwich panel
        %ilam=17;
        %ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=1;
        
    case 98 % Analysis: Cantilevered Unswept Wing y-dir span: Ribs & Spars at Midplane , Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 10;
        Mg = 10;

        panID=98;
        panXY = [ 0.0  0.0;  % Unswept rectangular wing y-dir
            50.0  0.0;
            50.0 110.0;
            0.0  110.0];   
        
        rsLayout=1234;
        nRib=6;
        nSpar=3;
        sparWebLamID=0;
        sparCapLamID=0;
        sparCapXsect=1;

        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=20; % sandwich panel
        %ilam=17;
        %ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=1;

    case 99 % Analysis: Cantilevered Unswept Wing y-dir span: Ribs at Midplane, Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 10;
        
        panID=99;
        panXY = [ 0.0  0.0;  % Unswept rectangular wing y-dir
            50.0  0.0;
            50.0 110.0;
            0.0  110.0];
        rsLayout=1234;
        nRib=0;
        nSpar=3;
        sparWebLamID=0;
        sparCapLamID=0;
        sparCapXsect=1;
       
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alumn
        %ilam=20; % sandwich panel
        %ilam=17;
        %ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=1;
        
                % Design data
        drespType = 'WEIGHT';
        % ... Thickness, Linear in y
%         tDeg=1;
%         XL=[0.02  -1.0  0.0  0.0];
%         XU=[1.50   1.0  0.0  0.0];

        % ... Thickness, Quadratic in y    
        %F = [1.0, y, y^2, x, x.*y ,y^2*x ,x^2 ,x^2*y ,y^2*x^2];
        tDeg=2;
        XL=[ 0.005 -1.00  -1.00  0.00  0.00  0.00   0.00   0.00  0.00];
        XU=[ 1.50   1.00   1.00  0.00  0.00  0.00   0.00   0.00  0.00];

        % .. Thickness, quadratic in x & y
%         tDeg=2;
%         XL=[ 0.005 -1.00  -1.00  -1.00  -1.00  -1.00   -1.00   -1.00  -1.00];
%         XU=[ 0.75   1.00   1.00   1.00   1.00   1.00    1.00    1.00   1.00];

        X0=0.5*(XL+XU);
        thk_coeff=X0;
        
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs
        
    case 100 % Analysis: Cantilevered VSP WING : Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 10;
        
       %panID=4;
        panID=100;
        dgfile = 'C:\Users\hensomc\Documents\VSP\OpenVSP-3.2.1-win32\work\learning\wing_numSect1_numU2.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);        
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        %ilam=0; % single layer alum
        ilam=20; % sandwich panel
        %ilam=21; % wing skewed laminate
        p0=-0.025; % uniform pressure
        plot_flag=1;
        
    case 101 % Analysis: Cantilevered Trapezoidal Wing : Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 10;
        Mg = 12;
        
        panID=12;
        panXY = [ 0.0  0.0;  % Wing Trapezoid
            70.0  0.0;
            70.0 10.0;
            0.0  25.0]; 
        
        bcType='cfff'; %'cfff'; %'ssss'; 'cccc';
        ilam=13;
        %ilam=20; % sandwich panel
        ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=0;
        
    case 102 % Analysis: Cantilevered Unswept Wing x-dir span: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 8;
        Mg = 10;
        
        panID=102;
        panXY = [ 0.0  0.0;  % Unswept rectangular wing x-dir
            110.0  0.0;
            110.0 50.0;
            0.0   50.0];
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=20; % sandwich panel
        %ilam=17;
        %ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=1;        
        
    % DELETE THIS CASE
    case 103 % Analysis: Unswept VSP Wing y-dir span: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        numMode = 15;
        ssCalc=1;
        pDeg = 6;
%         Mg = 8;
%         Mg = 16;
        Mg=8;
        
        panID=103;
        dgfile = 'unswept_wing_w_spars_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);               
        
        rsLayout=1234;
        nRib = 0;
        nSpar = 3;
        sparWebLamID=0;
        sparCapLamID=0;
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=20; % sandwich panel
        ilam=2;  % 8 ply quasi
        ilam=10; %32 ply quasi
        ilam=0; % alum
        %ilam=17;
        %ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=1;        
        
    case 104 % Analysis: SWEPT VSP Wing y-dir span, ribs & spars: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 10;
        
        panID=104;
        dgfile = 'rhsWing_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
                
        rsLayout=1234;
        nRib = 0;
        nSpar = 5;
        sparCapXsect=1;
        sparWebXsect=1;
        sparWebLamID=0;
        sparCapLamID=0;
        
        %bayMatName='nomex';
            
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alum
        %ilam=17;
        ilam=22;
        
%          p0=-0.025; % uniform pressure
        
        % Tip Torque Case
        ptLds = [   0.   0.    1.0;
                    0.   0.   -1.0 ];
        ptLoc = [  64.354    108.;
                   72.354    108.];

        % Tip In-plane Shear
%         ptLds = [   10.   0.    0.0;
%                     10.   0.    0.0 ];
%         ptLoc = [  64.354    108.;
%                    72.354    108.];
               
        plot_flag=0;
        
    case 105 % Analysis: Unswept VSP Wing y-dir span: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        numMode = 15;
        ssCalc=1;
        pDeg = 8;
        Mg=10;
        
        panID=105;
  
        %dgfile = 'vsp_unswept_rect_foil_rhs_DegenGeom.m';
        %dgfile = 'vsp_unswept_rect_foil_rhs_h0_05_DegenGeom.m';
        %dgfile = 'vsp_unswept_rect_foil_rhs_h0_6_DegenGeom.m';
        dgfile = 'vsp_unswept_rect_foil_rhs_h3_0_DegenGeom.m';                                                       
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);        
        
        rsLayout=1234;
        nRib = 0;
        nSpar = 3;
        sparCapXsect=0;
        sparWebXsect=1;
        sparWebLamID=0;
        sparCapLamID=0;
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        %ilam=20; % sandwich panel
        %ilam=2;  % 8 ply quasi
        %ilam=10; %32 ply quasi
        ilam=0; % alum
        %ilam=17;
        ilam=22;
        
        % Tip Point Load Case
%         ptLds = [   0.   0.   -1.0];
%         ptLoc = [  18.   108.];

        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [   9.   108.;
%                    27.   108.];
%         
        % Uniform Pressure Case
         p0=-0.025; % uniform pressure

        plot_flag=0;
        
    case 106 % Analysis: Cantilevered Trapezoidal Wing y-dir Liu Ref: Linear tapered thickness - Free Vib.
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 10;
        Mg = 12;
        
        panID=106;
        panXY = [ 0.0  0.0;  % Wing Trapezoid Liu Ref
            72.0   0.0;
            146.852  192.0;
            110.852 192.0];        
       
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0;
        %ilam=20; % sandwich panel
        
        tDeg=1; % Linear in x
        thk_coeff=[0.945 -0.855 0.0 0.0]; 

        % Uniform Pressure Case
        p0=-0.1; % uniform pressure

        % Tip torque Case
%         ptLds = [   0.   0.   1.0;
%                     0.   0.  -1.0];
%         ptLoc = [   110.852 192.0;
%                     146.852 192.0];
                
        plot_flag=1;

    case 107 % Analysis: Unswept Solid Fill VSP Wing y-dir span: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        numMode = 15;
        ssCalc=1;
        pDeg = 8;
        Mg=10;
        
        panID=105;
        dgfile = 'vsp_unswept_rect_foil_rhs_h3_0_DegenGeom.m';                                                       
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);        
       
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alum
        ilam=22; % vlam T0=0, T1=-30
        ilam=23; % vlam T0=0, T1=-30
        bayMatName='nomex';
        % Tip Point Load Case
%         ptLds = [   0.   0.   -1.0];
%         ptLoc = [  18.   108.];

        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [   0.0   108.;
%                    36.0   108.];
%         
        % Uniform Pressure Case
        p0=-0.5; % uniform pressure

        plot_flag=0;
        
    case 108 % Analysis: Unswept SOLID Fill VSP Wing w/ 3SPARS y-dir span: Ritz to FEM comparison - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        numMode = 15;
        ssCalc=1;
        pDeg = 6;
        Mg=10;

        panID=105;
        dgfile = 'vsp_unswept_rect_foil_rhs_h3_0_DegenGeom.m';                                                       
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);        
        
        rsLayout=1234;
        nRib = 3;
        nSpar = 3;
        sparWebLamID=0;
        sparCapLamID=0;
        sparCapXsect=1;
        sparWebXsect=0;
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alum
        bayMatName='nomex';
        % Tip Point Load Case
%         ptLds = [   0.   0.   -1.0];
%         ptLoc = [  18.   108.];

        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [   0.0   108.;
%                    36.0   108.];
%         
        % Used beam width = 6"
               
        % Uniform Pressure Case
        p0=-0.025; % uniform pressure

        plot_flag=1;
 
    case 109 % Analysis: F16 VSP COARSE Wing y-dir span, spars:  - Cantilevered with uniform pressure
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 12;
        
        panID=109;
        dgfile = 'vsp_f16_rhs_wing_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);   

        rsLayout=1234;
        nRib = 0;
        nSpar = 7;
        sparWebLamID=0;
        sparCapLamID=0;
        sparCapXsect=1;
        sparWebXsect=1;
        
        %bayMatName='syncore';
            
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alum
        %ilam=17;
        ilam=22;
        p0=-0.025; % uniform pressure
        plot_flag=1;

    case 110 % VSP LIU SOLID WING
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        numMode = 10;
        ssCalc=1;
        pDeg = 10;
        smearKM=1;
        Mg=12;
        
        panID=110;
        dgfile = 'vsp_rhsLiuSolidWing_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);  
 
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        %ilam=0; % vlam T0=0, T1=-30
        ilam=23; % vlam T0=0, T1=-30
        bayMatName='nomex';

        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [ 110.851   192.;
%                   146.851   192.];
%         
        % Uniform Pressure Case
       p0=-0.025; % uniform pressure

        plot_flag=1;
        
    case 111 % VSP LIU BUILT-UP WING, 4 SPARS, 10 RIBS
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        numMode = 15;
        ssCalc=1;
        pDeg = 6;
        Mg=12;
        smearKM=1;

        panID=110;
        dgfile = 'vsp_rhsLiuSolidWing_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);  

        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % vlam T0=0, T1=-30
  
        rsLayout=1234;
        nRib = 0;
        nSpar = 4;
        sparWebLamID=0;
        sparCapLamID=0;
        sparCapXsect=1;
        sparWebXsect=1;
        
        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [ 110.851   192.;
%                   146.851   192.];
        
        % Uniform Pressure Case
        p0=-0.025; % uniform pressure

        plot_flag=0;
        
    case 112 % VSP LIU BUILT-UP WING w/ SMEARED UNI LAM, 4 SPARS, 10 RIBS
        iDesOpt=0;
        isoltype=3;
        femSol=1;
        numMode = 10;
        ssCalc=1;
        pDeg = 6;
        Mg=10;
        smearKM=1;

        panID=110;
        dgfile = 'vsp_rhsLiuSolidWing_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);  

        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=23; % 4 ply smeared quasi
  
        rsLayout=1234;
        nRib = 0;
        nSpar = 4;
        sparWebLamID=0;
        sparCapLamID=0;
        sparCapXsect=1;
        sparWebXsect=1;
        
        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [ 110.851   192.;
%                   146.851   192.];
        
        % Uniform Pressure Case
        p0=-0.50; % uniform pressure

        plot_flag=1;

    case 113 % AIAA VSP LIU BUILT-UP WING w/ SMEARED 4 PLY VLAM, 4 SPARS, 10 RIBS
        
        iDesOpt=1;
        isoltype=3;
        femSol=1;
        numMode = 15;
        ssCalc=1;  % 1 - static cases
        pDeg = 4;  % 10 def
        Mg=30;   % 12 def
        smearKM=1;

        panID=110;
        dgfile = 'vsp_rhsLiuSolidWing_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);  

        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=22; % 4 ply smeared VLAM
  
%         tDeg=1; % Linear in x
%         thk_coeff=[0.945 -0.855 0.0 0.0]; 
        
        rsLayout=1234;
        nRib = 0;
        nSpar = 4;
        sparWebLamID=0; % alum
        sparCapLamID=0; % alum
        sparCapXsect=1;
        sparWebXsect=1;

        
        
        % Tip Torque Case
              
%         ptLds = 10000*[   0.   0.    1.0;
%                   0.   0.   -1.0 ];
%         ptLoc = [ 110.851   192.;
%                   146.851   192.]; 
%           Load acting at spar location
%           ptLoc = [ 118.051   192.;
%                   139.651   192.];
        
        % Tip Unit Bending Load - Flexural Axis
%         ptLds = -1000*[   0.   0.    1.0;
%                     0.   0.    1.0 ];
%         ptLoc = [ 110.851   192.;
%                   146.851   192.];     
%          ptLoc = [ 118.051   192.;
%                   139.651   192.];
% 
        % Tip Axial Load - Y dir
%         ptLds =  15000*[  0.   1.0    0.0;
%                           0.   1.0    0.0;
%                           0.   1.0    0.0;
%                           0.   1.0    0.0];
%         ptLoc = [ 118.051   192.;
%                   125.251   192.;
%                   132.451   192.;
%                   139.651   192.];
              
        % Uniform Pressure Case
        p0=-1; % uniform pressure

        plot_flag=1;
        
%         % Design data
%         drespType = 'EIGN';
%         XL=[-45 -45];
%         XU=[45 45];
%         X0=[20 45];
%         Xstart=cell(length(1));
%         Xstart{1} = X0;
%         
%         dvType=[1 1]; % Layer fiber angle DVs
        
        % Design data
        drespType = 'WEIGHT';
        % ... DV: Constant layer thickness only
%         XL=[0.005 0.005 0.005 0.005];
%         XU=[0.250 0.250 0.250 0.250];
%         dvType=[0];
        
        % ... DV: Constant layer Thickness + Theta
%         XL=[0.005 0.005 0.005 0.005 -90 -90];
%         XU=[0.250 0.250 0.250 0.250  90  90];   
%         dvType=[0 1];

        % ... DV: Constant thickness - code mod needed to support this
%         tDeg=0;
%         XL=[0.02   0.0  0.0  0.0];
%         XU=[1.50   0.0  0.0  0.0];
%         dvType=[2];
        
        % ... DV: Variable total thickness, linear in y
%         tDeg=1;
%         XL=[0.02  -1.0  0.0  0.0];
%         XU=[1.50   1.0  0.0  0.0];
%         dvType=[2];

        % ... DV: Variable total thickness, quadratic in y
        %F = [1.0, y, y^2, x, x.*y ,y^2*x ,x^2 ,x^2*y ,y^2*x^2];
%         tDeg=2;
%         XL=[ 0.02  -1.00  -1.00  0.00  0.00  0.00   0.00   0.00  0.00];
%         XU=[ 1.50   1.00   1.00  0.00  0.00  0.00   0.00   0.00  0.00];
%         dvType=[2];

        %... DV: Variable total thickness, quadratic in x, y
        %F = [1.0, y, y^2, x, x.*y ,y^2*x ,x^2 ,x^2*y ,y^2*x^2];
%         tDeg=2;
%         XL=[ 0.02  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00];
%         XU=[ 1.50   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00];
%         dvType=[2];

        %... DV: Theta + Variable total thickness, quadratic in x, y
        %F = [1.0, y, y^2, x, x.*y ,y^2*x ,x^2 ,x^2*y ,y^2*x^2];
        tDeg=2;
        XL=[ -90 -90  0.02  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00  -1.00 ];
        XU=[  90  90  1.50   1.00   1.00   1.00   1.00   1.00   1.00   1.00   1.00];
        dvType=[1 2];
                
        % ... Start points
        X0=0.5*(XU+XL);
        % ... Near optimium
        %X0=[0.0478   -0.0801    0.0   0.0];
        %X0=[0.0478   -0.0801    0.0442         0         0         0         0         0         0];

        if tDeg>0
            %thk_coeff=X0;
            thk_coeff=X0(1,3:end);
        end
        Xstart=cell(length(1));
        Xstart{1} = X0;

    case 114 % AIAA Sci Tech 2017 Uniform Trapezoidal Plate y-dir, Linear tapered thickness
        iDesOpt=1;
        isoltype=3;
        femSol=1;
        numMode = 10;
        ssCalc=1;
        pDeg = 5;
        Mg = 12;

        panID=114;
        panXY = [ 0.0  0.0;  % AIAA SciTech 2017 Uniform Trapezoid
            72.0   0.0;
            54.492  192.0;
            18.492 192.0];

        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0;
        %ilam=20; % sandwich panel
        
%         tDeg=1; % Linear in y
%         thk_coeff=[0.945 -0.855 0.0 0.0]; 
       
        % Uniform Pressure Case
        p0=-0.05; % uniform pressure

        % Tip torque Case
%         ptLds = [   0.   0.   1.0;
%                     0.   0.  -1.0];
%         ptLoc = [   18.492 192.0;
%                     54.492 192.0];       
        plot_flag=0;
        
        % Design data
        drespType = 'WEIGHT';
        % ... Thickness, Linear in y
%         tDeg=1;
%         XL=[0.02  -1.0  0.0  0.0];
%         XU=[1.50   1.0  0.0  0.0];

        % ... Thickness, Quadratic in y    
        %F = [1.0, y, y^2, x, x.*y ,y^2*x ,x^2 ,x^2*y ,y^2*x^2];
        tDeg=2;
        XL=[ 0.005 -1.00  -1.00  0.00  0.00  0.00   0.00   0.00  0.00];
        XU=[ 1.50   1.00   1.00  0.00  0.00  0.00   0.00   0.00  0.00];

        % .. Thickness, quadratic in x & y
%         tDeg=2;
%         XL=[ 0.005 -1.00  -1.00  -1.00  -1.00  -1.00   -1.00   -1.00  -1.00];
%         XU=[ 0.75   1.00   1.00   1.00   1.00   1.00    1.00    1.00   1.00];

        X0=0.5*(XL+XU);
        thk_coeff=X0;
        
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs
 
    case 115 % AIAA SciTech 2017 Uniform Trapezoidal Core Filled Wing, Rect Foil
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        numMode = 10;
        ssCalc=1;
        pDeg = 10;
        Mg=12;

        panID=115;
        dgfile = 'C:\Users\hensomc\Documents\VSP\MyWork\rhsUniformTrapezoidRectFoil_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);

        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alum
        ilam=23; % QI
        bayMatName='syncore';
        % Tip Point Load Case
%         ptLds = [   0.   0.   -1.0];
%         ptLoc = [  18.   108.];

        % Tip Torque Case
        ptLds = [   0.   0.    1.0;
                    0.   0.   -1.0 ];
        ptLoc = [   18.492 192.0;
                    54.492 192.0];

%         
        % Uniform Pressure Case
%         p0=-0.5; % uniform pressure

        plot_flag=0;
        
      case 120 % ESAVE Wing
        iDesOpt=0;
        isoltype=1;
        femSol=1;
        numMode = 10;
        ssCalc=1;
        pDeg = 5;
        Mg=12;

        panID=120;
        dgfile = 'vsp_esav_wing_rhs_DegenGeom.m';
        [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
        
        bcType='fcff'; %'cfff'; %'ssss'; 'cccc';
        ilam=0; % alum
        %ilam=23; % QI
        %ilam=22; % 4 ply smeared VLAM
        % Substructure
        rsLayout=1234;
        nRib = 0;
        nSpar = 19;
        sparWebLamID=0; % alum
        sparCapLamID=0; % alum
        sparCapXsect=1;
        sparWebXsect=1;
        
        % Tip Point Load Case
%         ptLds = [   0.   0.   -1.0];
%         ptLoc = [  18.   108.];

        % Tip Torque Case
%         ptLds = [   0.   0.    1.0;
%                     0.   0.   -1.0 ];
%         ptLoc = [   18.492 192.0;
%                     54.492 192.0];

%         
        % Uniform Pressure Case
%         p0=-0.5; % uniform pressure

        plot_flag=1;      
        
     % -------------------------
     % DESIGN OPTIMIZATION CASES
     % -------------------------
     case 200 % Design Opt, Cantilevered plate, uniform pressure, Minimize Max Strain, Weight constraint = 2, bi-material laminate
        iDesOpt=1;
        isoltype=3;
        femSol=1;
        ssCalc=1;
        pDeg = 5;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'cfff'; %'ssss'; 'cccc';
        ilam=7; % Bi-material alum, c/ep;
        p0=-0.025; % uniform pressure
        plot_flag=1;

        % Design data
        drespType = 'STRAIN';
        XL=[0.005 0.005];
        XU=[0.250 0.250];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;
        dvType=[0]; % Layer thickness DVs

        
     case 201 % Design Opt, cfff plate, Max vib Freq. Obj, Max Weight Constraint = 2, bi-material laminate
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=0;
        pDeg = 5;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'cfff'; %'ssss'; 'cccc';
        ilam=7; % Bi-material alum, c/ep;
        plot_flag=1;

        % Design data
        drespType = 'EIGN';
        XL=[0.005 0.005];
        XU=[0.250 0.250];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;
        dvType=[1]; % Layer thick DVs

    
     case 202 % Design Opt, ssss plate, Min Weight Obj, Strain Constraint, bi-material laminate
        iDesOpt=1;
        isoltype=2;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        p0=-0.025; % uniform pressure
        bcType='ssss'; %'ssss'; 'cccc'; 'cfff';
        ilam=7; % Bi-material alum, c/ep;
        plot_flag=1;
        
        % Design data
        drespType = 'WEIGHT';
        XL=[0.005 0.005];
        XU=[0.250 0.250];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;
        dvType=[0]; % Layer thick DVs
        
    case 203 % Design Opt, ssss plate, Min Weight Obj, Constraints: [Buckling, Strain], bi-material laminate
        iDesOpt=1;
        isoltype=2;
        femSol=1;
        ssCalc=1;
        pDeg = 5;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='ssss'; %'ssss'; 'cccc'; 'cfff';
        ilam=7;  % 13, 14, 15, 2, 9
        %Nx=-100; Ny=   0; Nxy=0; % Nx buckling case
        Nx=  0; Ny=-100; Nxy=0; % Ny buckling case
        ilam=7; % Bi-material alum, c/ep;
        plot_flag=1;        
 
        % Design data
        drespType = 'WEIGHT';
        XL=[0.005 0.005];
        XU=[0.250 0.250];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;
        dvType=[0]; % Layer thick DVs
        
    case 204 % Design Opt, ssss plate, Min Weight Obj, Constraints: [Strain, LINEAR varTHK],  - Al laminate
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 5;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'ssss'; 'cccc'; 'cfff';
        p0=-0.025; % uniform pressure
        ilam=0; % 1 ply alum;
        plot_flag=1;        

        % Design data
        drespType = 'WEIGHT';
        tDeg=1; % Linear in x
        XL=[0.02   0.0 -0.15  0.0];
        XU=[0.15   0.0  0.15  0.0];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs


    case 205 % Design Opt, ssss plate, Min Weight Obj, Constraints: [Strain, QUADRATIC varTHK],  - Al laminate
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'ssss'; 'cccc'; 'cfff';
        p0=-0.025; % uniform pressure
        ilam=0; % 1 ply alum;
        plot_flag=1;        

        % Design data
        drespType = 'WEIGHT';
        tDeg=2; % Quadratic in x
        XL=[ 0.005  0.00  0.00  -0.05  0.00  0.00  -0.05   0.00  0.00];
        XU=[ 0.15   0.00  0.00   0.15  0.00  0.00   0.15   0.00  0.00];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        %dvType=[2 2 2 2 2 2 2 2 2]; % Poly coeffs
        dvType=[2]; % Poly coeffs
        
    case 206 % Design Opt, ssss plate, Min Weight & Max Freq Obj, Constraints: [Strain],  - bi-material laminate
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'ssss'; 'cccc'; 'cfff';
        p0=-0.025; % uniform pressure
        ilam=7;  % bi-metal
        plot_flag=1;        

        % Design data
        drespType = 'WEIGHT_EIGN';
        XL=[0.005 0.005];
        XU=[0.250 0.250];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[0]; % thk

    case 207 % Design Opt, ssss plate, Min Weight & Max Freq Obj, Linear Var Thk., Constraints: [Strain],  - Al laminate
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 20;

        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'ssss'; 'cccc'; 'cfff';
        p0=-0.025; % uniform pressure
        ilam=0;  % bi-metal
        plot_flag=1;        

        % Design data
        drespType = 'WEIGHT_EIGN';
        tDeg=1;
        XL=[0.02   0.0 -0.15  0.0];
        XU=[0.15   0.0  0.15  0.0];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs


    case 208 % Design Opt, ssss plate, Min Weight & Max Freq Obj, Quadratic thickness, Constraints: [Strain],  - bi-material laminate
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 20;
 
        panID=4;
        panXY = [ 0.0  0.0;  % Square
            10.0  0.0;
            10.0 10.0;
            0.0 10.0];

        bcType='cfff'; %'ssss'; 'cccc'; 'cfff';
        p0=-0.025; % uniform pressure
        ilam=0;  % bi-metal
        plot_flag=1;        

        % Design data
        drespType = 'WEIGHT_EIGN';
        tDeg=2; % Quadratic in x
        XL=[ 0.005  0.00  0.00  -0.05  0.00  0.00  -0.05   0.00  0.00];
        XU=[ 0.15   0.00  0.00   0.15  0.00  0.00   0.15   0.00  0.00];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs

      
    case 209 % Design Opt, Min Weight Obj, Constraints: [Strain, LINEAR varTHK], Geom: Cantilevered Trapezoidal Wing y-dir Liu Ref: Linear tapered thickness
        iDesOpt=1;
        isoltype=1;
        femSol=1;
        ssCalc=1;
        pDeg = 6;
        Mg = 20;

        panID=106;       
        panXY = [ 0.0  0.0;  % Wing Trapezoid Liu Ref
            72.0   0.0;
            146.852  192.0;
            110.852 192.0];
        
        bcType='fcff'; %'ssss'; 'cccc'; 'cfff';
        ilam=0;  %
        plot_flag=1;        

        % Design data
        drespType = 'WEIGHT';
        tDeg=1;
        %thk_coeff=[0.945 -0.855 0.0 0.0]; orig thickness profile
        XL=[0.02  -1.5  0.0  0.0];
        XU=[0.945  0.0 0.0 0.0];
        %XU=[0.15   0.0  0.15  0.0];
        X0=0.5*(XU+XL);
        Xstart=cell(length(1));
        Xstart{1} = X0;      
        thk_coeff=X0;
        dvType=[2]; % Poly coeffs

        
        % Uniform Pressure Case
        p0=-0.025; % uniform pressur

        % Tip torque Case
%         ptLds = [   0.   0.   1.0;
%                     0.   0.  -1.0];
%         ptLoc = [   110.852 192.0;
%                     146.852 192.0];
        plot_flag=1;
        
    otherwise
        disp('Invalid jobID');   
end

% Zero polythk coeffs
if(tDeg==0)
    thk_coeff=[0 0 0 0];
end;

% Panel geom and layout
%[panXY, degenGeom]=getPanGeom(panID, degenGeomFileName, plot_flag);

% Revise Mg if needed for Spars & Ribs
Mg = check_Mg(nSpar,Mg);

% Set return arguments
solDATA = setSolData( ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc, paramType );
%lamDATA = getLamInputs(ilam); eliminate?
lamDATA = setLamData(ilam, tDeg, thk_coeff);  
ldsDATA = setLdsData( Nx, Ny, Nxy, ptLds, ptLoc, p0 );

[ribDATA,sparDATA,panels]=genRibSpar(rsLayout, nRib, nSpar, panXY);

% Combine structures into wingDATA
wingDATA.panXY = panXY;

wingDATA.spar = sparDATA;
wingDATA.sparCapXsect = sparCapXsect;
wingDATA.sparWebXsect = sparWebXsect;
wingDATA.sparWebLamID = sparWebLamID;
wingDATA.sparCapLamID = sparCapLamID;

wingDATA.degenGeom = degenGeom;

wingDATA.rib = ribDATA;

wingDATA.panels = panels;

wingDATA.bayMatName = bayMatName;

desDATA=setDesData(XL,XU,X0,Xstart,Rallow,drespType,dvType);   


function [solDATA]=setSolDefaults()
%% Set solver defaults

iDesOpt=0;       % 0 = analysis only, 1 = design optimization, 2= parameter study
isoltype=1;      % 1=vib,    2=buckl,    3=linear static
ipoltype=2;      % 1=power,  2=legendre, 3=power-all numerical
intMethod=2;     % 1=simpson rule, 2=gauss quadrature
RitzSE=0;        % 1=use Ritz super-element
symSol=0;        % 1=check using symbolic solution
femSol=0;        % 1=create a Nastran FE model for solution
ssCalc=0;        % 1=recover streass and strain results
%vspDegenGeom=0;  % 1=import data from OpenVSP DgenGeom file
iThkFun=0;       % 1=use polynomial thickness function

bcMethod=2;      % 1=spring, 2=null space

pDeg=8;         % initial polynomial degree (use ODD value for isoltype = 1)

smearKM=0;       % 0=discrete ply&mass stiffness integration, 1 = smeared at z=z0

Mg=20;           % number of integration pts - psi dir, assume Ng=Mg (use even value to get nodes at midspan)

numMode=15;      % Number of modes to recover

paramType=0;

plot_flag=0;     % 0 = none, 1=modes/info, 2=animate modes/debug, 3=all details

% BC inputs
bcType='ssss';      % 4 edge BC's of panel, 'c'=clamped, 's'=simple, 'f'=free
bcSpring = 1.0e10;  % use: 1e10 = vibration; 1e13 = buckling|statics

% Set solution paramaters
solDATA = setSolData( ipoltype,pDeg,Mg,smearKM,plot_flag,isoltype,iDesOpt,femSol,symSol,intMethod,bcMethod,bcType,RitzSE,bcSpring,iThkFun,numMode,ssCalc, paramType );


function [lamDATA] = getLamInputs(ilam)
%% Get laminate input data

            %ilam = 0  % 1-ply isotropic alum plate
            %ilam = 1  % 2-ply isotropic layup
            %ilam = 2  % 8-ply Quasi layup
            %ilam = 3  % 2-ply curvilinear layup <+-45|0>
            %ilam = 4  % 8-ply antisymmetric curvilinear layup <+-45|0>
            %ilam = 5  % 4-ply antisymmetric [+45/-45/+45/-45]
            %ilam = 6  % 2-ply bi-metallic aluminum-steel laminate [0 0]
            %ilam = 7  % 2-ply aluminum-carbon_epoxy [0 0]
            %ilam = 8  % 8-ply Uni layup
            %ilam = 9  % 8-ply Angle ply +/-45 layup
            %ilam =10  % 32-ply Quasi layup
            %ilam =11  %  2-ply +/-45 layup
            %ilam =12  %  2-ply +/-30 layup
            %ilam =13  %  8-ply curvilinear layup:
            %[<-45|45>,<45|-45>,<-45|45>,<45|-45>,  <45|-45>,<-45|45>,<45|-45>,<-45|45>]
            %ilam =14  % Gurdal 2008-Case I laminate: 8 ply curvilinear, %[<0|45>,<0|-45>,<0|45>,<0|-45>,  <0|-45>,<0|45>,<0|-45>,<0|45>]
            %ilam =15: % Gurdal 2008-Case II laminate: 8 ply curvilinear, %[<90|45>,<90|-45>,<90|45>,<90|-45>,  <90|-45>,<90|45>,<90|-45>,<90|45>]
            %ilam =16: 32-ply Angle Ply layup
            %ilam =17: 8 ply curvilinear, %[<45|0>,<-45|0>,<45|0>,<-45|0>,  <-45|0>,<45|0>,<-45|0>,<45|0>]

% Thickness coefficients [c0, c1, c2, c3, c4]; t=c0*t0 +c1*x +c2*y + c3*xy
%tDeg=0;
% tDeg=1;
%thk_coeff=[0.1   0.0   0.0   0.0]; % constant thickness
% thk_coeff=[0.15   0.0   -0.1/2   0.0]; % linear variation x-dir
%thk_coeff=[0.1   0.1     0.1   0.1]; % bi-linear variation xy-dir

%Laminate data structure
%[lamDATA] = setLamData(ilam, tDeg, thk_coeff);  


% function [panXY, degenGeom]=getPanGeom(panID,dgfile,plot_flag)
% %% Get panel input data
% % Panel geometry
% 
% % Default vaule
% degenGeom=[];
% 
% switch panID
%     
%     case 0 % DegenGeom file
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
%         
%     case 1 
%         panXY = [-1.0 -1.0;  % Normalized Square
%             1.0 -1.0;
%             1.0  1.0;
%             -1.0  1.0];
%     case 2 
%         panXY = [-5.0 -5.0;  % Centered Square
%             5.0 -5.0;
%             5.0  5.0;
%             -5.0  5.0];
%     case 3
%         panXY = [ 0.0  0.0;  % Strip
%             10.0  0.0;
%             10.0  1.0;
%             0.0  1.0];
%     case 4
%         panXY = [ 0.0  0.0;  % Square
%             10.0  0.0;
%             10.0 10.0;
%             0.0 10.0];
%     case 5
%         panXY = [ 0.0  0.0;  % Large Square
%             20.0  0.0;
%             20.0 20.0;
%             0.0 20.0];
%     case 6
%         panXY = [ 0.0  0.0;  % Wing Trapezoid
%             70.0  0.0;
%             70.0 10.0;
%             0.0  25.0];
%     case 7
%         panXY = [ 0.0  0.0;  % Square - Houmat ref
%             11.8  0.0;
%             11.8 11.8;
%             0.0 11.8];
%     case 8
%         panXY = [ 0.0  0.0;  % x-bias Rectangle
%             20.0  0.0;
%             20.0 10.0;
%             0.0 10.0];
%     case 9
%         panXY = [ 1.0  2.0;  % Quadrilateral
%             11.0  2.0;
%             9.0 12.0;
%             1.0 10.0];
%     case 10
%         panXY = [-4.0 -4.0;  % Quadrilateral, same as above, shifted
%             6.0 -4.0;
%             4.0  6.0;
%             -4.0  4.0];
%     case 11
%         panXY = [ -20.00          0.00  % Convex Quad
%             40.00          8.00
%             20.00         25.00
%             -35.00         14.00];
%     case 12       
%         panXY = [ 0.0  0.0;  % Wing Trapezoid
%             70.0  0.0;
%             70.0 10.0;
%             0.0  25.0]; 
%    case 13       
%         panXY = [ 0.0  0.0;  % Dr Wang test plate
%             30.0  0.0;
%             30.0 20.0;
%             0.0  20.0];    
%     case 14
%         panXY = [ 0.0  0.0;  % y-bias Rectangle
%             10.0  0.0;
%             10.0 20.0;
%             0.0 20.0];
%     case 98
%         panXY = [ 0.0  0.0;  % Unswept rectangular wing y-dir
%             50.0  0.0;
%             50.0 110.0;
%             0.0  110.0];   
%     case 99
%         panXY = [ 0.0  0.0;  % Unswept rectangular wing y-dir
%             50.0  0.0;
%             50.0 110.0;
%             0.0  110.0];
% 
%     case 100
%         dgfile = 'C:\Users\hensomc\Documents\VSP\OpenVSP-3.2.1-win32\work\learning\wing_numSect1_numU2.m';
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
%         
%     case 102
%         panXY = [ 0.0  0.0;  % Unswept rectangular wing x-dir
%             110.0  0.0;
%             110.0 50.0;
%             0.0   50.0];
%     case 103
%         dgfile = 'unswept_wing_w_spars_DegenGeom.m';
%         %dgfile = 'vsp_unswept_rhswing_9pts_DegenGeom.m';
%         
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
%         
%     case 104
%         dgfile = 'rhsWing_DegenGeom.m';
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
%         
%     case 105
%         %dgfile = 'vsp_unswept_rect_foil_rhs_DegenGeom.m';
%         %dgfile = 'vsp_unswept_rect_foil_rhs_h0_05_DegenGeom.m';
%         %dgfile = 'vsp_unswept_rect_foil_rhs_h0_6_DegenGeom.m';
%         dgfile = 'vsp_unswept_rect_foil_rhs_h3_0_DegenGeom.m';                                                       
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
%         
%     case 106
%         panXY = [ 0.0  0.0;  % Wing Trapezoid Liu Ref
%             72.0   0.0;
%             146.852  192.0;
%             110.852 192.0];
% 
%     case 109
%          dgfile = 'vsp_f16_rhs_wing_DegenGeom.m';
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);   
%         
%     case 110
%          dgfile = 'vsp_rhsLiuSolidWing_DegenGeom.m';
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);  
% 
%     case 114
%         panXY = [ 0.0  0.0;  % AIAA SciTech 2017 Uniform Trapezoid
%             72.0   0.0;
%             54.492  192.0;
%             18.492 192.0];
%         
%     case 115  % AIAA Uniform Trapezoidal Core Filled Wing, Rect Foil
%         dgfile = 'C:\Users\hensomc\Documents\VSP\MyWork\rhsUniformTrapezoidRectFoil_DegenGeom.m';
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
% 
%     case 120  % ESAVE Wing
%         dgfile = 'vsp_esav_wing_rhs_DegenGeom.m';
%         [degenGeom,panXY] = importDegenGeom(dgfile,plot_flag);
%         
%     otherwise
%         disp('Invalid PANEL ID')
% end



function [ldsDATA]=setLdsDefaults()
%% Get load input data

% Running loads
Nx=0; Ny=0; Nxy=0;

% Point Loads
ptLds=0; ptLoc=0;
% ptLds = [   0.   0.   -0.5;
%             0.   0.   -0.5 ];
% ptLoc = [  10.   0.;
%            10.  10.];

% ptLds = [   0.   0.   -0.5 ];
% ptLoc = [  10.   5.];

% Distributed pressure loads
p0=0.0;
%p0 = -0.025;

       
% ptLoc = [  10.   0.;  % For a unit strip
%            10.   1.];
            
[ldsDATA] = setLdsData(Nx,Ny,Nxy, ptLds, ptLoc, p0); % Nx,Ny,Nxy,...


% function [solDATA,lamDATA,panXY,ldsDATA] = readProblem(probID,solDATA,lamDATA,panXY,ldsDATA)
% %% Get problem inputs
% switch probID
%     case 1  % Vib frequency convergence: QI laminate, SSSS
%         solDATA.soltype=1;
%         solDATA.bcType='ssss';
%         lamDATA.ilam=2;
%         paramSweep = 1; 
% end
% 
% 
% function [solDATA,lamDATA,panXY,ldsDATA]=setDefaults(solDATA,lamDATA,panXY,ldsDATA)
% %% Set default problem inputs
% 
% % get solver inputs
% solDATA=getSolInputs();
% 
% % get laminate data
% lamDATA=getLamInputs();
% 
% % get panel geometry
% panXY=getPanGeom();
% 
% % get loads data
% ldsDATA=getLdsInputs();
% 

function [XL,XU,X0,Xstart,drespType,Rallow,dvType]=setDesDefaults()
%% Get Design Optimization inputs

    % DESIGN DATA
    % Side constraints and initial values
    
    Rallow = 12; % allowable fiber placement turning radius
    
    XL=[];
    XU=[];
    X0=[];
    Xstart=cell(length(1));
    drespType='';
    dvType=[];
    
    %drespType = 'EIGN';   % response type
    %drespType = 'WEIGHT';
    %drespType = 'BUCKL';
    %drespType = 'STRAIN';

    % Legendre thickness polynomial:
    %   t2 = y.^2;
    %   t3 = t2.*(3.0./2.0);
    %   t4 = t3-1.0./2.0;
    %   t5 = x.^2;
    %   t6 = t5.*(3.0./2.0);
    %   t7 = t6-1.0./2.0;
    %
    %   F = [1.0, y,         t4, x, x.*y,          t4.*x,         t7,          t7.*y,                    t4.*t7];
    %   F = [1.0, y, 1.5y*y-0.5, x,  x*y, (1.5y*y-0.5)*x, 1.5x^2-0.5, (1.5x^2-0.5)*y, (1.5y*y-0.5)*(1.5x*x-0.5)];

    % Constant thickness - original
    %XL=[0.005 0.005];
    %XU=[0.250 0.250];
    
    % Fiber angles only: [<T0|T1>]
%     XL=[-45 -45];
%     XU=[ 45  45];    
    %XL=[  0 -90];
    %XU=[ 90  90];    

%     % Constant thickness + fiber angles
%     XL=[0.005 0.005 -45 -45];
%     XU=[0.250 0.250  45  45];
    
    % Linear in x: tDeg=1
    %XL=[0.02   0.0 -0.15  0.0];
    %XU=[0.15   0.0  0.15  0.0];

    % Bilinear in x & y: tDeg=1
%     XL=[-0.05  -0.05 -0.05   -0.05];
%     XU=[0.15    0.15  0.15    0.15];

    % Linear in x:
%     XL=[-1.00   0.00  0.00  -1.00  0.00  0.00   0.00   0.00  0.00];
%     XU=[ 1.00   0.00  0.00   1.00  0.00  0.00   0.00   0.00  0.00];

    % Bi-Linear in x: tDeg=2
%     XL=[-1.00  -1.00  0.00  -1.00 -1.00  0.00   0.00   0.00  0.00];
%     XU=[ 1.00   1.00  0.00   1.00  1.00  0.00   0.00   0.00  0.00];

      % Quadratic in x: Neg thk
%     XL=[-1.00   0.00  0.00  -1.00  0.00  0.00  -1.00   0.00  0.00];
%     XU=[ 1.00   0.00  0.00   1.00  0.00  0.00   1.00   0.00  0.00];
      
      % Quadratic in x: Pos thk
%     XL=[-0.05   0.00  0.00  -0.05  0.00  0.00  -0.05   0.00  0.00];
%     XU=[ 0.15   0.00  0.00   0.15  0.00  0.00   0.15   0.00  0.00];

      % Quadratic in x & y:
%     XL=[-0.05  -0.05 -0.05  -0.05  0.00  0.00  -0.05   0.00  0.00];
%     XU=[ 0.15   0.15  0.15   0.15  0.00  0.00   0.15   0.00  0.00];   

%     XL=[-1.00  -1.00 -1.00  -1.00  0.00  0.00  -1.00   0.00  0.00];
%     XU=[ 1.00   1.00  1.00   1.00  0.00  0.00   1.00   0.00  0.00];   

      % Complete Quadratic in x & y
%     XL=[-0.05  -0.05 -0.05  -0.05 -0.05 -0.05  -0.05  -0.05 -0.05];
%     XU=[0.15    0.15  0.15   0.15  0.15  0.15   0.15   0.15  0.15];   
%     X0=[ 0.10   0.00  0.00   0.00  0.00  0.00   0.00   0.00  0.00];    
    
    % Use for constant thickness
    %X0=XL;
    %X0=0.5*(XU+XL);
    %X0=[45 0];
    
    %Xstart=cell(length(1));  % Best choices:vertices & center of feasible region (how can we quickly find them? 
    %Xstart{1} = X0;
%    Xstart{1} = [ 45   0];
%     Xstart{2} = [  5 -85];
%     Xstart{3} = [ 85 -85];
%     Xstart{4} = [ 85  85];
%     Xstart{5} = [  5  85];

    % Use for variable thickness - will need to return lamDATA, or modify
%      ilam = 0;
%     tDeg=2;
%      thk_coeff=X0;
%      [lamDATA] = setLamData(ilam,tDeg,thk_coeff);

function [ ribDATA, sparDATA, panels] = genRibSpar( rsLayout, nRib, nSpar, panXY )
%genRibSpar generates rib and spar definitions based on rsLayout parameter
%   rsLayout = (12, 34, 1234) = rib spar layout
%            = 12 spars parallel to side 12 (trailing edge)
%            = 34 spars parallel to side 34 (leading edge)
%            = 1234 spars "fanned" out at equal chord fractions between
%              LE and TE

% Set return arguments with default values
ribDATA=[];
ribDATA.XY=[];

sparDATA=[];
sparDATA.XY=[];
panels=[];
if (rsLayout <= 0)
    return
end

switch rsLayout
    case 12 % Spars parallel to T.E.
        
    case 34 % Spars parallel to L.E.
        
    case 1234 % Spars fanned out equally fro L.E. to T.E.
        %sparDATA2.XY=cell(2,2,3)
        rootChordX = panXY(2,1)-panXY(1,1);
        tipChordX  = panXY(3,1)-panXY(4,1);
        rootChordY = panXY(2,2)-panXY(1,2);
        tipChordY  = panXY(3,2)-panXY(4,2);
        
        sparSpacRoot = rootChordX/(nSpar+1);
        sparSpacTip  = tipChordX/(nSpar+1);
        for i=1:nSpar          
            xspar1 = panXY(1,1)+ i/(nSpar+1) * rootChordX;
            xspar2 = panXY(4,1)+ i/(nSpar+1) * tipChordX;
            
            yspar1 = panXY(1,2)+ i/(nSpar+1) * rootChordY;
            yspar2 = panXY(4,2)+ i/(nSpar+1) * tipChordY;

            %sparDATA.XY=[sparDATA.XY; xspar1 yspar1; xspar2,yspar2];
            sparDATA.XY{i}=[xspar1 yspar1; xspar2,yspar2];
        end
        
        % Ribs
        for i=1:nRib
            ribDATA.XY{1}=[panXY(1,1) panXY(1,2); panXY(2,1) panXY(2,2)];
            ribDATA.XY{2}=[panXY(4,1) panXY(4,2); panXY(3,1) panXY(3,2)];
            ribDATA.XY{3}=[panXY(4,1) panXY(4,2); panXY(3,1) panXY(3,2)];
        end
    otherwise
        disp('Invalid rib and spar layout parameter = '), rsLayout
end

% Geometry

% Material and section properties
sparDATA.sectType(1)=1;
sparDATA.matType(1)=1;

% Extract Skin Panels from Rib/Spar Layout
[panels]=getSkinPanels(nRib, nSpar, panXY, ribDATA, sparDATA);



%% Get Skin Panels
function [skinPanels]=getSkinPanels(nRib, nSpar, panXY, ribDATA, sparDATA)

skinPanels=[];

% Return if no sub-structure defined
if nSpar==0 && nRib == 0
    return
end

spar=sparDATA;
rib=ribDATA;

% Concatenate spar endpts with panel points
nKeyPts=nSpar*2+length(panXY); % assumes no ribs
nRegion=(nSpar+1)*(nRib+1);
keypts = linspace(1,nKeyPts,nKeyPts);

% Number of panel nodes along Wing x & y
if nSpar==0
    nnx=1;
else
    nnx=nSpar+2;
end

if nRib==0
    nny=1;
else
    nny=nRib+2;
end

% 
sparXY=[];
for i=1:nSpar
    sparXY=[sparXY;cell2mat(spar.XY(i))];
end
XYData = [keypts' [panXY;sparXY]]

% Build matrix of panel regions
Region=[];
if nSpar==0
    Region = [1 keypts nnx nny];
else
    snode1=5;
    snode2=6;
    for i=1:nRegion
        if(i==1) %1st region
            Region = [Region; i XYData(1,1,1) XYData(5,1,1) XYData(6,1,1) XYData(4,1,1) nnx nny];
        elseif (i == nRegion)
            Region = [Region; i XYData(nKeyPts-1,1,1) XYData(2,1,1) XYData(3,1,1) XYData(nKeyPts,1,1) nnx nny];
        else
            Region = [Region; i XYData(snode1,1,1) XYData(snode1+2,1,1) XYData(snode2+2,1,1) XYData(snode2,1,1) nnx nny];
            snode1=snode1+2;
            snode2=snode2+2;
        end
    end
end
% Region
% for i=1:nRegion
%     disp('Panel: ');i
%     for j=1:4
%         disp('cornerPts='); XYData(Region(i,j+1),2), XYData(Region(i,j+1),3)
%     end
% end

panels.Regions=Region;
panels.XY=XYData;
    



