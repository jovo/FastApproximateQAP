function [f0, g ] = fungrad(x,A,B,C,alpha)
%function [f0, g ] = fungrad(x,A,B)
% gradient (all matrix level computations)
%
% Louis J. Podrazik circa 1996
%
% IDA Center for Computing Sciences
%  (c) 1996, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%
    
    % In http://arxiv.org/pdf/1112.5507v5.pdf it looks as though the
    % objective function is tr(A P B^T P^T) but here it suggests it is
    % tr(B^T Q A^T P^T) = tr(P A Q^T B)
    [m,n]=size(A);
    [P,Q]=unstack(x,m,n);
    %f0=sum(sum(P*A*Q'.*B));
    f0=fun(x,A,B,C,alpha);
    %Note sum(sum((P*A*Q').*B))=sum(sum(B*Q*A').*P)
    %  so grad wrt to P is  B*Q*A'
    g=reshape(2*(1-alpha)*B*Q*A' + alpha*C,m*m,1);
    %Note sum(sum((P*A*Q').*B))=sum(sum(B'*P*A).*Q)
    %  so grad wrt to Q is  B'*P*A
    %g=-[g;reshape(B'*P*A,n*n,1)];
    g=[g;reshape(2*(1-alpha)*B'*P*A + alpha*C,n*n,1)];

