%% generate .mat file from txt elegans connectomes
clear, clc

rootDir='~/Research/projects/primary/FastApproximateQAP/';

%% get Achem and permutations
Achem = dlmread(sprintf([rootDir '/data/graphm/data/Achem.txt']));
[nvertices foo]=size(Achem);

permed_Achem=zeros(nvertices,nvertices,100);
for idx=1:100
    permed_Achem(:,:,idx) = dlmread(sprintf([rootDir '/data/graphm/data/Achem_permd/Achem_%d.txt'],idx));
end

%% get Agap and permutations
Agap = dlmread(sprintf([rootDir '/data/graphm/data/Agap.txt']));
[nvertices foo]=size(Achem);

permed_Agap=zeros(nvertices,nvertices,100);
for idx=1:100
    permed_Agap(:,:,idx) = dlmread(sprintf([rootDir '/data/graphm/data/Agap_permd/Agap_%d.txt'],idx));
end

%% get actual permuations

perms=zeros(nvertices,100);
iv = 1:nvertices;
for idx=1:100
    the_perm = dlmread(sprintf([rootDir, '/data/graphm/data/perm/perm_%d.txt'],idx));
    the_perm = the_perm(:); 
    perms(the_perm,idx) = iv;
end

%%
save([rootDir, 'data/elegans_connectomes2.mat'])

