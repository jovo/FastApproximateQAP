function [s,p,iters,bests,t,n]=run_qaplib_problems(problems,niter,r)
%function [s,p,iters,bests,t,n]=run_qaplib_problems(problems,niter,r)
% Run the list of QAPLIB test problems given in the cell array problems
% using at most niter iteration of the Frank-Wolf method (FW). On output
% s, a length(problems) long array of final scores.
% p, a length(problems) long cell array permutations that achieve the final
%     final scores given in the array s.
% iters, a vector given the number of iterations of FW that were used.
% bests, if known the current best known answer for problem i will be
%        returned in bests(i).
% t,    is a vector of elapsed times for problems.
% n,    gives the problems sizes for problems.
% John M. Conroy, June 2012
%
% IDA Center for Computing Sciences
%  (c) 2012, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%
s=zeros(length(problems),1); iters=s; bests=s; t=s; n=s;
p=cell(size(s));
for i=1:length(problems)
    tic;
    [A,B,~,bests(i)]=qap_read(problems{i},'data/qaplib');
    if (nargin==3)&&(~isempty(r))
        [s(i),p{i},~,iters(i)]=sfw(A,B,niter,-1);
    else
        [s(i),p{i},~,iters(i)]=sfw(A,B,niter);    
    end
    t(i)=toc;
    n(i)=length(p{i});
end

