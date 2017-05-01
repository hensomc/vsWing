function getDispErrFEM(panXY, XYData1, solDATA, lamDATA, nasDisp ,ritzDisp ,plot_flag)

% Performs comparison of Ritz and FEM displacement vectors

isoltype = solDATA.isoltype;


%% Error
%eigError = abs((nasEigVals - EgN(1:10)./EgN(1:10)))*100;

% Table of max/min displacements
% format long;
% disp('=============================================')
% disp('Max/Min displacement comparison: FEM and Ritz');
% disp('=============================================')
% index=min([numEig,length(nasEigVals),length(nasEigVals)]);
% T1=table( nasEigVals(1:index), ritzEigVals(1:index), ...
%          'VariableNames',{'nasEigVals','ritzEigVals'});

ritzMax = [max(max(ritzDisp(:,1))); max(max(ritzDisp(:,2))); max(max(ritzDisp(:,3))); max(max(ritzDisp(:,4))); max(max(ritzDisp(:,5)))];
ritzMin = [min(min(ritzDisp(:,1))); min(min(ritzDisp(:,2))); min(min(ritzDisp(:,3))); min(min(ritzDisp(:,4))); min(min(ritzDisp(:,5)))];
nasMax  = [max(max(nasDisp(:,1))); max(max(nasDisp(:,2))); max(max(nasDisp(:,3))); max(max(nasDisp(:,4))); max(max(nasDisp(:,5)))];
nasMin  = [min(min(nasDisp(:,1))); min(min(nasDisp(:,2))); min(min(nasDisp(:,3))); min(min(nasDisp(:,4))); min(min(nasDisp(:,5)))];
rowhdr  = ['Disp_u'; 'Disp_v'; 'Disp_w'; 'Phi_x '; 'Phi_y '];

table_title='Max/Min Displacement Summary';
disp(table_title);
T1=table (rowhdr, ritzMax, nasMax, ritzMin, nasMin);
format short
disp(T1)

titles={'ritzMax', 'nasMax', 'ritzMin', 'nasMin'};
a=[ritzMax ritzMin nasMax nasMin];
%table2word(titles,a,'vlamWordOutputTables.doc',table_title);

% Section cut @ y=b/2, values along x-axis
midpt1=0.5*[ panXY(4,1)+panXY(1,1)  panXY(4,2)+panXY(1,2)]
midpt2=0.5*[ panXY(2,1)+panXY(3,1)  panXY(2,2)+panXY(3,2)]
xydata2pt=[1  midpt1;
           2  midpt2]
sectn1=Find_Line_Node([1 2], xydata2pt, XYData1)
xvals=XYData1(sectn1, 2)

% Section cut @ x=a/2, values along y-axis
midpt1=0.5*[ panXY(1,1)+panXY(2,1)  panXY(1,2)+panXY(2,2)];
midpt2=0.5*[ panXY(3,1)+panXY(4,1)  panXY(3,2)+panXY(4,2)];
xydata2pt=[1  midpt1;
           2  midpt2];
sectn2=Find_Line_Node([1 2], xydata2pt, XYData1);

% LE Section cut @ x=0
xydata2pt=[1  panXY(1,1) panXY(1,2);
           2  panXY(4,1) panXY(4,2)];
sectnLE=Find_Line_Node([1 2], xydata2pt, XYData1);

% TE Section cut @ x=a
xydata2pt=[1  panXY(2,1) panXY(2,2);
           2  panXY(3,1) panXY(3,2)];
sectnTE=Find_Line_Node([1 2], xydata2pt, XYData1);


% Get title strings
[panTitle,polTitle,lamTitle,propTitle,lamPropTitle] = getPanTitles(panXY, solDATA, lamDATA);

% Graph section1 displacements

% w-displacements
figure
plot( xvals, ritzDisp(sectn1,3), 'r:o','linewidth',2), hold on
plot( xvals, nasDisp(sectn1,3), 'b-.s','linewidth',2), hold off

title( {'';'\bfMidspan Displacement w as Function of (x)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation x (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Midspan Displacement w as Function of (x)');
T=table (xvals, ritzDisp(sectn1,3), nasDisp(sectn1,3));
disp(T)

% v-displacements
figure
plot( xvals, ritzDisp(sectn1,2), 'r:o','linewidth',2), hold on
plot( xvals, nasDisp(sectn1,2), 'b-.s','linewidth',2), hold off

title( {'';'\bfMidspan Displacement v as Function of (x)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation x (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Midspan Displacement v as Function of (x)');
T=table (xvals, ritzDisp(sectn1,2), nasDisp(sectn1,2));
disp(T)


% Graph section2 displacements
% w-displacements
xvals=XYData1(sectn2, 3)

figure
plot( xvals, ritzDisp(sectn2,3), 'r:o','linewidth',2), hold on
plot( xvals, nasDisp(sectn2,3), 'b-.s','linewidth',2), hold off

title( {'';'\bfMidspan Displacement w as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Midspan Displacement w as Function of (y)');
T=table (xvals, ritzDisp(sectn2,3), nasDisp(sectn2,3));
disp(T)

% v-displacements
figure
plot( xvals, ritzDisp(sectn2,2), 'r:o','linewidth',2), hold on
plot( xvals, nasDisp(sectn2,2), 'b-.s','linewidth',2), hold off

title( {'';'\bfMidspan Displacement v as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('Ritz','FEM','location','best');

disp('Midspan Displacement v as Function of (y)');
T=table (xvals, ritzDisp(sectn2,2), nasDisp(sectn2,2));
disp(T)

% LE & TE w-displacements
xvals=XYData1(sectnLE, 3);
figure
plot( xvals, ritzDisp(sectnLE,3), 'r:o','linewidth',2), hold on
plot( xvals, nasDisp(sectnLE,3), 'b-.s','linewidth',2), hold on

disp('LE Displacement w as Function of (y)');
T=table (xvals, ritzDisp(sectnLE,3), nasDisp(sectnLE,3));
disp(T)

xvals=XYData1(sectnLE, 3);
plot( xvals, ritzDisp(sectnTE,3), 'g:d','linewidth',2), hold on
plot( xvals, nasDisp(sectnTE,3), 'b-.x','linewidth',2), hold off

title( {'';'\bfLE & TE Displacement w as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('LE-Ritz','LE-FEM','TE-Ritz','TE-FEM','location','best');

disp('TE Displacement w as Function of (y)');
T=table (xvals, ritzDisp(sectnLE,3), nasDisp(sectnLE,3));
disp(T)

% Combined Midspan u,v,w displacments
xvals=XYData1(sectn2, 3);
figure
plot( xvals, 10*ritzDisp(sectn2,1), 'r:o','linewidth',2), hold on
plot( xvals, 10*nasDisp(sectn2,1), 'r-.s','linewidth',2), hold on

plot( xvals, 10*ritzDisp(sectn2,2), 'g:o','linewidth',2), hold on
plot( xvals, 10*nasDisp(sectn2,2), 'g-.s','linewidth',2), hold on

plot( xvals, ritzDisp(sectn2,3), 'b:o','linewidth',2), hold on
plot( xvals, nasDisp(sectn2,3), 'b-.s','linewidth',2), hold on
title( {'';'\bfMidspan Displacements u,v,w as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('10u Ritz','10u FEM','10v Ritz','10v FEM','w Ritz','w FEM','location','best');


% Combined LE u,v,w displacments
xvals=XYData1(sectnLE, 3);
figure
plot( xvals, 100*ritzDisp(sectnLE,1), 'r:o','linewidth',2), hold on
plot( xvals, 100*nasDisp(sectnLE,1), 'r-.s','linewidth',2), hold on

plot( xvals, 100*ritzDisp(sectnLE,2), 'g:o','linewidth',2), hold on
plot( xvals, 100*nasDisp(sectnLE,2), 'g-.s','linewidth',2), hold on

plot( xvals, ritzDisp(sectnLE,3), 'b:o','linewidth',2), hold on
plot( xvals, nasDisp(sectnLE,3), 'b-.s','linewidth',2), hold on
title( {'';'\bfLE Displacements u,v,w as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('100u Ritz','100u FEM','100v Ritz','100v FEM','w Ritz','w FEM','location','best');

% Combined TE u,v,w displacments
xvals=XYData1(sectnTE, 3);
figure
plot( xvals, 100*ritzDisp(sectnTE,1), 'r:o','linewidth',2), hold on
plot( xvals, 100*nasDisp(sectnTE,1), 'r-.s','linewidth',2), hold on

plot( xvals, 100*ritzDisp(sectnTE,2), 'g:o','linewidth',2), hold on
plot( xvals, 100*nasDisp(sectnLE,2), 'g-.s','linewidth',2), hold on

plot( xvals, ritzDisp(sectnTE,3), 'b:o','linewidth',2), hold on
plot( xvals, nasDisp(sectnTE,3), 'b-.s','linewidth',2), hold on

title( {'';'\bfTE Displacements u,v,w as Function of (y)';panTitle;lamTitle;'';''},'fontsize',14 )
xlabel('\bfLocation y (in.)','fontsize',12)
ylabel('\bf Displacement (in.)','fontsize',12)
legend('100u Ritz','100u FEM','100v Ritz','100v FEM','w Ritz','w FEM','location','best');