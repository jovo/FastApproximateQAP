function P=perm2mat(p)
%function P=perm2mat(p)

% Louis J. Podrazik circa 1996
%
% IDA Center for Computing Sciences
%  (c) 1996, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

n=max(size(p));
P=zeros(n,n);
for i=1:n
   P(i,p(i))=1;
end;
