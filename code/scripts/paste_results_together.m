%% paste together results from eric harley

%
% err(i,j,1) == error from algorithm i on run j w/ chem
% err(i,j,2) == error from algorithm i on run j w/ gap

clear, clc
rootDir='~/Research/projects/primary/FastApproximateQAP';
rootDir0=rootDir;
load([rootDir, '/data/results/elegans_results_graphm_eric.mat']); % load FAQ results

allErrors=1-err;
allTimes=time;

%% add FAQ Chem
rootDir=rootDir0;
load([rootDir, '/data/results/elegans_results_FAQChem_eric.mat']); % load FAQ results

allErrors(1,:,1)=errors(1:999)/279;
allTimes(1,:,1)=times(1:999);

%% add FAQ Gap
rootDir=rootDir0;
load([rootDir, '/data/results/elegans_results_FAQGap_eric.mat']); % load FAQ results

allErrors(1,:,2)=errors(1:999)/279;
allTimes(1,:,2)=times(1:999);

%% save

save('../../data/results/elegans_results_all.mat','allErrors','allTimes')