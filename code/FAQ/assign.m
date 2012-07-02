function [p,w,x]=assign(A,munk)
%function [p,w,x]=assign(A,munk)
% The optional second argument if set to false will use 
% the maxassign_linprog() the
% interior point method based on matlab's optimization toolbox.
% munkres(), an alternating path method, is used by default.
%
% John M. Conroy 1996-2010
% conroy@super.org
% IDA Center for Computing Sciences
%  (c) 1996-2010, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%
if nargin==1,
    munk=true;
end


if munk==false
    [p,w,x]=maxassign_linprog(A'); p=p';
else
    [p,w]=lapjv(-A',0.01); w=-w;
     %   [p,w]=munkres(-A'); w=-w;
    x=perm2mat(p)';
end
%  x=reshape(x,size(A));
%
%[q,T]=hungarian(-A);
% if length(unique(p))~=length(p)
%     [p,c]=munkres(-A');
%     x=perm2mat(p);
% end
% if norm(pm-p')~=0,
%    % disp('Non-unique')
%     p=pm;
% end
%x=perm2mat(p)';
%


% OLD CODE used for 1999 experiments and calls an external binary.
% n=size(A,2);
% mmA=max(max(abs(A)));
% A=1.0e4/mmA*A;
% A=round(max(max(A))-A);
% %[t,p]=assct(A);
% fid=fopen('/tmp/A.mat','w');
% fprintf(fid,'%d \n',n);
% for j=1:n
% for i=1:n
%    fprintf(fid,'%d \n',A(i,j));
% end
% end
% !/home/conroy/matlab/qap/assign
% load -ascii /tmp/p.mat
