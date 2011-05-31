function [A,s]=stoch(A,dim)
% function A=stoch(A,dim)
% Normalize A so that it is row stochastic if
% one input argument is given or along the
% dim dimension 1 for column stochastic and
% 2 for row stochastic.
% Currently only 2D arrays are supported.
%
% John M. Conroy circa 1996
%
% IDA Center for Computing Sciences
%  (c) 1996, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

[m,n]=size(A);
if (nargin==1)|(dim==2)
    s=full(sum(A,2));
    if any(s)==0
        s=max(s,realmin);
        warning('Zero sum found!');
    end
    A=spdiags(1./s,0,m,m)*A;
else
    s=full(sum(A,1))';
    if any(s)==0
        s=max(s,realmin);
        warning('Zero sum found!');
    end
    A=A*spdiags(1./s,0,n,n);
end
