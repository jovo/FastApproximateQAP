function [f0] = fun(x,A,B,C,alpha)
% function [f0] = fun(x,A,B)
%--------------------
% Louis J. Podrazik circa 1996
%
% IDA Center for Computing Sciences
%  (c) 1996, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

[m,n]=size(A);
[P,Q]=unstack(x,m,n);
f0 = 2*(1-alpha)*trace(P*A*Q'*B) + alpha * trace(C'*P);  % Labelled graph matching equation
%f0=-f0; % MINIMIZATION
