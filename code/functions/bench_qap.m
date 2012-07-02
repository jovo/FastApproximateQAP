function [times,errors,iters]=bench_qap(data,Nperms)
% The Fast Approximate Quadratic Assignment Problem (FAQAP) algorithm is
% run comparing some input square matrix with a randomly permuted version.
% This comparison is run either on an input graph or a randomly generated
% Erdos-Renyi graph.
% 
% INPUT:
%   data: could be an adjacency matrix or a vector
%       -- if data is a vector, it is interpreted as a vector of digraph sizes.  
%           Each digraph is ER with p=log(n)/n.
%       -- if data is a square matrix, is it interpreted as the input, and
%           we run it both letting data be a weighted graph, and binarizing
%       -- else, terminate: this code requires square matrices
%   Nperms: number of random permutations to run, only if input is a square
%           matrix
% 
% OUTPUT:
%   times: run time per iteration
%   errors: the number of edges in the permuted graph that disagree
%   iterations: number of iterations of FAQAP performed is given in iters.
% 
% 
% If the data is a graph, on return the three arguments times, errors, and iters will be as before
% except they will have two columns each one for the first experiment and
% one for the experiment with the entries of A set to 0 and 1.
% 
% John M. Conroy
% September 9, 2011
% IDA Center for Computing Sciences
%  (c) 1996-2010, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%
% modified by Joshua Vogelstein, 2012

s = RandStream('mt19937ar','Seed',0);
[m,n]=size(data);
if (n>1)&&(m==n)
    if varargin==1, Nperms=10; end
    times=zeros(Nperms,2);
    errors=times; iters=errors;
    for i=1:Nperms
        q=randperm(s,n);
        disp(i)
        for j=1:2
            if j==1,
                A=data;
            else
                A=spones(data);
            end
            B=A;
            B(q,q)=A;
            tic;
            [~,myq,~,iters(i,j)]=sfw(A,-B,30);
            times(i,j)=toc;
            errors(i,j)=sum(myq~=q);
        end
    end
else
    times=zeros(size(data));
    errors=times;
    iters=errors;
    parfor i=1:length(data)
        
        display(['numher of vectices = ', num2str(i)]) 
        % Generate a random digraph is Pr(i,j)=log(n)/n;
        p=log(data(i))/data(i);
        A=spones(sprand(data(i),data(i),p));
        q=randperm(s,data(i));
        B=A;
        B(q,q)=B;
        tic;
        [~,myq,~,iters(i)]=sfw(A,-B,30);
        times(i)=toc;        
        errors(i)=sum(myq~=q);
    end
end
