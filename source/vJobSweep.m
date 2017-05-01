function vJobSweep()
%vJobSweep Sweeps over a collection of vsWing jobs
%   Detailed explanation goes here

%dgfileList={'vsp_rhsWing_DegenGeom.m'}


% Path to VSP design files - uncomment addpath call below
%addpath('C:\Users\hensomc\Documents\MATLAB\vs_lam\vspDesignSpaceFiles\');
%dgfileList={'LWD_SW_0_AR_0_TR_0_DegenGeom.m'};

dgfileList={
 'LWD_SW_2_AR_0_TR_0_DegenGeom.m',
 'LWD_SW_2_AR_0_TR_1_DegenGeom.m',
 'LWD_SW_2_AR_0_TR_2_DegenGeom.m',
 'LWD_SW_2_AR_1_TR_0_DegenGeom.m',
 'LWD_SW_2_AR_1_TR_1_DegenGeom.m',
 'LWD_SW_2_AR_1_TR_2_DegenGeom.m',
 'LWD_SW_2_AR_2_TR_0_DegenGeom.m',
 'LWD_SW_2_AR_2_TR_1_DegenGeom.m',
 'LWD_SW_2_AR_2_TR_2_DegenGeom.m',
 'LWD_SW_0_AR_0_TR_0_DegenGeom.m',
 'LWD_SW_0_AR_0_TR_1_DegenGeom.m';
 'LWD_SW_0_AR_0_TR_2_DegenGeom.m',
 'LWD_SW_0_AR_1_TR_0_DegenGeom.m',
 'LWD_SW_0_AR_1_TR_1_DegenGeom.m',
 'LWD_SW_0_AR_1_TR_2_DegenGeom.m',
 'LWD_SW_0_AR_2_TR_0_DegenGeom.m',
 'LWD_SW_0_AR_2_TR_1_DegenGeom.m',
 'LWD_SW_0_AR_2_TR_2_DegenGeom.m',
 'LWD_SW_1_AR_0_TR_0_DegenGeom.m',
 'LWD_SW_1_AR_0_TR_1_DegenGeom.m',
 'LWD_SW_1_AR_0_TR_2_DegenGeom.m',
 'LWD_SW_1_AR_1_TR_0_DegenGeom.m',
 'LWD_SW_1_AR_1_TR_1_DegenGeom.m',
 'LWD_SW_1_AR_1_TR_2_DegenGeom.m',
 'LWD_SW_1_AR_2_TR_0_DegenGeom.m',
 'LWD_SW_1_AR_2_TR_1_DegenGeom.m',
 'LWD_SW_1_AR_2_TR_2_DegenGeom.m'
};

% Summary figure
figure;

%ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')

ha=tight_subplot(3,9,[.01 .03],[.1 .01],[.01 .01]);
for i=1:length(dgfileList)
    [wingDATA]=vMain(0, dgfileList{i});
    
    % draw quad in x-y and corner node IDs
    XY4=wingDATA.panXY
    xy=[XY4;XY4(1,:)];
%     ax(i)=subplot(3,9,i);
    axes(ha(i));
    plot(xy(:,1),xy(:,2),'k','linewidth',2)
	
end
set(ha(:),'XTickLabel',''); set(ha(:),'YTickLabel','');
linkaxes(ha,'xy');
axis off
