clear, clc
rootDir='~/Research/projects/primary/FastApproximateQAP';
rootDir0=rootDir;
load([rootDir, '/data/results/elegans_results_faq_ericComp.mat']); % load FAQ results
rootDir=rootDir0;
load([rootDir, '/data/results/elegans_results_graphm_eric.mat']); % load FAQ results
rootDir=rootDir0;
savestuff=1;

%% scatter plots
% Have to add FAQ runs.
%
% err(i,j,1) == error from algorithm i on run j w/ chem
% err(i,j,2) == error from algorithm i on run j w/ gap
%
% time(i,j,1) == time of algorithm i on run j w/ chem
% time(i,j,2) == time of algorithm i on run j w/ gap
%
% perm(i,j,1,:) == the permutation produced by algorithm i on run j w/ chem
% perm(i,j,2,:) == the permutation produced by algorithm i on run j w/ gap

% algorithm indexes
% FAQ_idx  = 1;
% I_idx    = 2;
% rand_idx = 3;
% RANK_idx = 4;
% U_idx    = 5;
% QCV_idx  = 6;
% PATH_idx = 7;

% Machine specs:
%       Processor Name: Intel Core i7
%       Processor Speed: 2.2 GHz
%       Number of Processors: 1
%       Total Number of Cores: 4
%       L2 Cache (per Core): 256 KB
%       L3 Cache: 6 MB
%       Memory: 8 GB

markers{1}='s';
markers{2}='h';
markers{3}='+';
markers{4}='x';
markers{5}='o';
markers{6}='d';
markers{7}='^';

colors{1}='k';
colors{2}='y';
colors{3}='b';
colors{4}='m';
colors{5}='g';
colors{6}='b';
colors{7}='r';

algsKept=[1, 5:7];

markersizes=12*ones(7,1);

err(1,1:100,1)=1-chem.errors/279;
err(1,1:100,2)=1-gap.errors/279;

time(1,1:100,1)=chem.times;
time(1,1:100,2)=gap.times;

n = 279;



%% plot results
clc

left=0.1;
bottom1=0.55;
bottom2=0.1;
width=0.80;
height=0.3;

%% error
figure(1); clf
subplot(221) %'Position',[left bottom1 width height]);
boxplot(1-err(algsKept,1:100,1)'...
    ,'labels',{'FAQ','U','QCV','PATH'}...
    ,'labelorientation','horizontal'...
    ,'color','k','symbol','k+','plotstyle','compact'); %...
%     ,'positions',linspace(0,.1,4));
title('Chemical')
ylabel('error (%)')
set(gca,'YTick',0:.5:1)
ylim([0 1])

subplot(222) %'Position',[left bottom2 width height]);
boxplot(1-err(algsKept,1:100,2)'...
    ,'labels',{'FAQ','U','QCV','PATH'}...
    ,'labelorientation','horizontal'...
    ,'color','k','symbol','k+','plotstyle','compact');
ylim([0 1])
set(gca,'YTick',0:.5:1)
title('Electrical')

%% time
% figure(1); clf
subplot(223), cla %'Position',[left bottom1 width height]);
boxplot(time(algsKept,1:100,1)'...
    ,'labels',{'FAQ','U','QCV','PATH'}...
    ,'labelorientation','horizontal'...
    ,'color','k','symbol','k+','plotstyle','compact'); %...
ylabel('time (seconds)')

subplot(224), cla %'Position',[left bottom2 width height]);
boxplot(time(algsKept,1:100,2)'...
    ,'labels',{'FAQ','U','QCV','PATH'}...
    ,'labelorientation','horizontal'...
    ,'color','k','symbol','k+','plotstyle','compact');

%%

if savestuff==1
    wh=[3 1.3]*2;
    figName=[rootDir, '/figs/connectomes'];
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print('-dpdf',figName)
    print('-dpng',figName)
    saveas(gcf,figName)
end

%%
percentiles=[0.05, .25, .5, .75, .95];
chemFAQTimePercentiles=quantile(time(1,1:100,1),percentiles)
chemFAQTimeIQR=chemFAQTimePercentiles(4)-chemFAQTimePercentiles(2)
