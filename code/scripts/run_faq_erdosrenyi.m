%% make cubic figure for QAP paper
clear, close all

ns = repmat(100:100:1000,1,1);
[times,errors,iters]=bench_qap(ns);

p=polyfit(ns,times,3);

n=10000;
np=p(1)*n^3+p(2)*n^2+p(3)*n+p(4);
nyears=np/60/60/24/365;

savestuff=0;
fname= '../data/ErdosRenyi_results.mat';
if savestuff, save(fname), end