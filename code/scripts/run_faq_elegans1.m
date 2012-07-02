%% connectome times
clear, clc

load('../data/elegans_connectomes.mat')

[chem.times,chem.errors,chem.iters]=bench_qap(Achem);
[gap.times,gap.errors,gap.iters]=bench_qap(Agap);

[median(chem.times); mad(chem.times)]
[median(chem.errors); mad(chem.errors)]
[median(chem.iters); mad(chem.iters)]

[median(gap.times); mad(gap.times)]
[median(gap.errors); mad(gap.errors)]
[median(gap.iters); mad(gap.iters)]

save('../data/elegans_results.mat')
