%% connectome times
clear, clc

load('~/Research/data/EM/c_elegans/ConnOrdered_040903.mat')

[times,errors,iters]=bench_qap(A_init_t_ordered);

[median(times); mad(times)]
[median(errors); mad(errors)]
[median(iters); mad(iters)]


[times,errors,iters]=bench_qap(Ag_t_ordered);

[median(times); mad(times)]
[median(errors); mad(errors)]
[median(iters); mad(iters)]
