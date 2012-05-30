%% connectome times
clear, clc

load('../data/elegansGraph.mat')

[chem.times,chem.errors,chem.iters]=bench_qap(Achem);

[median(chem.times); mad(chem.times)]
[median(chem.errors); mad(chem.errors)]
[median(chem.iters); mad(chem.iters)]


[gap.times,gap.errors,gap.iters]=bench_qap(Agap);

[median(gap.times); mad(gap.times)]
[median(gap.errors); mad(gap.errors)]
[median(gap.iters); mad(gap.iters)]

load('../data/elegansGraph_results.mat')
