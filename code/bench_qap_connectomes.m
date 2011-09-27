%% connectome times

load('~/Research/data/EM/c_elegans/c_elegans_chemical_connectome.mat')

[times,errors,iters]=bench_qap(A_init_t_ordered);

[median(times); mad(times)]
[median(errors); mad(errors)]
[median(iters); mad(iters)]


[times,errors,iters]=bench_qap(Ac);

[median(times); mad(times)]
[median(errors); mad(errors)]
[median(iters); mad(iters)]
