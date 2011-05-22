function [d,q] = dsproj(x,g,m,n)
%function [d,q] = dsproj(x,g,m,n)
% minimize gradient times d subject to linear constraints
% Louis J. Podrazik circa 1996
% Modified by John M. Conroy, 1996-2010
%
% IDA Center for Computing Sciences
%  (c) 1996-2010, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

[P,Q]=unstack(x,m,n);
[gP,gQ]=unstack(g,m,n);


[q,wq,wQ]=assign(-gQ);
wP=wQ;

w=stack(wP,wQ,m,n);
d=w-x;
