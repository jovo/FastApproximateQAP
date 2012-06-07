function A=sink(A,n)
%function A=sink(A,n)
% Perform n iterations of Sinkhorn balancing.
%
% John M. Conroy circa 1996
% conroy@super.org
% IDA Center for Computing Sciences
%  (c) 1996-2010, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

for i=1:n
    A=stoch(A,1);
    A=stoch(A,2);
end
