% Attached MAT file with three matrices: err, time, perm.
% 
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
% 
% algorithm indexes
% FAQ_idx  = 1;
% I_idx    = 2;
% rand_idx = 3;
% RANK_idx = 4;
% U_idx    = 5;
% QCV_idx  = 6;
% PATH_idx = 7;

% ---------------------------------------------
clear;

%%

n = 279;
num_runs = 999;

% FAQ I rand RANK U QCV PATH
num_algs = 7;

chem_idx = 1;
gap_idx = 2;

FAQ_idx  = 1;
I_idx    = 2;
rand_idx = 3;
RANK_idx = 4;
U_idx    = 5;
QCV_idx  = 6;
PATH_idx = 7;

%%

load 'output/Achem_U_times.txt'
load 'output/Achem_I_times.txt'
load 'output/Achem_PATH_times.txt'
load 'output/Achem_QCV_times.txt'
load 'output/Achem_RANK_times.txt'
load 'output/Achem_rand_times.txt'

load 'output/Agap_U_times.txt'
load 'output/Agap_I_times.txt'
load 'output/Agap_PATH_times.txt'
load 'output/Agap_QCV_times.txt'
load 'output/Agap_RANK_times.txt'
load 'output/Agap_rand_times.txt'

[Y,I] = sort(Achem_U_times(:,1));

Achem_U_times = Achem_U_times(I,:);
Achem_I_times = Achem_I_times(I,:);
Achem_PATH_times = Achem_PATH_times(I,:);
Achem_QCV_times = Achem_QCV_times(I,:);
Achem_RANK_times = Achem_RANK_times(I,:);
Achem_rand_times = Achem_rand_times(I,:);

Agap_U_times = Agap_U_times(I,:);
Agap_I_times = Agap_I_times(I,:);
Agap_PATH_times = Agap_PATH_times(I,:);
Agap_QCV_times = Agap_QCV_times(I,:);
Agap_RANK_times = Agap_RANK_times(I,:);
Agap_rand_times = Agap_rand_times(I,:);


%%

err  = zeros(num_algs,num_runs,2);
perm = zeros(num_algs,num_runs,2,n);
time = zeros(num_algs,num_runs,2);

for run_idx = 1:num_runs


    iv = 1:n;
    rp = zeros(n,1);
    the_perm = dlmread(sprintf('data/perm/perm_%d.txt',run_idx));
    the_perm = the_perm(:);
    rp(the_perm) = iv;
    

    graphm_perms = dlmread(sprintf('output/exp_out_file_%d',run_idx),' ',5,0);

    p_I_chem = graphm_perms(:,1);
    p_rand_chem = graphm_perms(:,5);
    p_RANK_chem = graphm_perms(:,3);
    p_U_chem = graphm_perms(:,2);
    p_QCV_chem = graphm_perms(:,4);
    p_PATH_chem = graphm_perms(:,6);

    % estimated permutation
    % perm(FAQ_idx,run_idx,chem_idx) = p_FAQ_chem;
    perm(I_idx,run_idx,chem_idx,:) = p_I_chem;
    perm(rand_idx,run_idx,chem_idx,:) = p_rand_chem;
    perm(RANK_idx,run_idx,chem_idx,:) = p_RANK_chem;
    perm(U_idx,run_idx,chem_idx,:) = p_U_chem;
    perm(QCV_idx,run_idx,chem_idx,:) = p_QCV_chem;
    perm(PATH_idx,run_idx,chem_idx,:) = p_PATH_chem;

    % error (# incorrectly assigned vertices)
    % err(FAQ_idx,run_idx,chem_idx) = sum( rp == p_FAQ_chem ) / n;
    err(I_idx,run_idx,chem_idx) = sum( rp == p_I_chem ) / n;
    err(rand_idx,run_idx,chem_idx) = sum( rp == p_rand_chem ) / n;
    err(RANK_idx,run_idx,chem_idx) = sum( rp == p_RANK_chem ) / n;
    err(U_idx,run_idx,chem_idx) = sum( rp == p_U_chem ) / n;
    err(QCV_idx,run_idx,chem_idx) = sum( rp == p_QCV_chem ) / n;
    err(PATH_idx,run_idx,chem_idx) = sum( rp == p_PATH_chem ) / n;

    graphm_perms = dlmread(sprintf('output/exp_out_file_gap_%d',run_idx),' ',5,0);

    p_I_gap = graphm_perms(:,1);
    p_rand_gap = graphm_perms(:,5);
    p_RANK_gap = graphm_perms(:,3);
    p_U_gap = graphm_perms(:,2);
    p_QCV_gap = graphm_perms(:,4);
    p_PATH_gap = graphm_perms(:,6);

    % estimated permutation
    % perm(FAQ_idx,run_idx,gap_idx) = p_FAQ_gap;
    perm(I_idx,run_idx,gap_idx,:) = p_I_gap;
    perm(rand_idx,run_idx,gap_idx,:) = p_rand_gap;
    perm(RANK_idx,run_idx,gap_idx,:) = p_RANK_gap;
    perm(U_idx,run_idx,gap_idx,:) = p_U_gap;
    perm(QCV_idx,run_idx,gap_idx,:) = p_QCV_gap;
    perm(PATH_idx,run_idx,gap_idx,:) = p_PATH_gap;

    % error (# incorrectly assigned vertices)
    % err(FAQ_idx,run_idx,gap_idx) = sum( rp == p_FAQ_gap ) / n;
    err(I_idx,run_idx,gap_idx) = sum( rp == p_I_gap ) / n;
    err(rand_idx,run_idx,gap_idx) = sum( rp == p_rand_gap ) / n;
    err(RANK_idx,run_idx,gap_idx) = sum( rp == p_RANK_gap ) / n;
    err(U_idx,run_idx,gap_idx) = sum( rp == p_U_gap ) / n;
    err(QCV_idx,run_idx,gap_idx) = sum( rp == p_QCV_gap ) / n;
    err(PATH_idx,run_idx,gap_idx) = sum( rp == p_PATH_gap ) / n;

end

for run_idx = 1:num_runs
    % wall time
    % time(FAQ_idx,run_idx,chem_idx)  = Achem_FAQ_times(run_idx,2);
    time(I_idx,run_idx,chem_idx)    = Achem_I_times(run_idx,2);
    time(rand_idx,run_idx,chem_idx) = Achem_rand_times(run_idx,2);
    time(RANK_idx,run_idx,chem_idx) = Achem_RANK_times(run_idx,2);
    time(U_idx,run_idx,chem_idx)    = Achem_U_times(run_idx,2);
    time(QCV_idx,run_idx,chem_idx)  = Achem_QCV_times(run_idx,2);
    time(PATH_idx,run_idx,chem_idx) = Achem_PATH_times(run_idx,2);

    % time(FAQ_idx,run_idx,gap_idx)   = Agap_FAQ_times(run_idx,2);
    time(I_idx,run_idx,gap_idx)     = Agap_I_times(run_idx,2);
    time(rand_idx,run_idx,gap_idx)  = Agap_rand_times(run_idx,2);
    time(RANK_idx,run_idx,gap_idx)  = Agap_RANK_times(run_idx,2);
    time(U_idx,run_idx,gap_idx)     = Agap_U_times(run_idx,2);
    time(QCV_idx,run_idx,gap_idx)   = Agap_QCV_times(run_idx,2);
    time(PATH_idx,run_idx,gap_idx)  = Agap_PATH_times(run_idx,2);

end