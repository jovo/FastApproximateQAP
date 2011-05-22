function [p,w,x]=maxassign_linprog(C)
%function [p,w,x]=maxassign_linprog(C)
% Find the maximum assignment of C
% calls matlab's linprog routine
% where p and w are such that
% w=0; for i=1:n, w=w+C(i,p(i)); end
% and x is the interior point solution given
% by linprog, which when reshaped to an n by n
% matrix should be close to a permutation matrix
% John M. Conroy circa 1996
%
% IDA Center for Computing Sciences
%  (c) 1996-2010, Institute for Defense Analyses, 4850 Mark Center Drive, Alexandria, Virginia, 22311-1882; 703-845-2500.
%
%     This material may be reproduced by or for the U.S. Government pursuant to the copyright license under the clauses at DFARS 252.227-7013 and 252.227-7014.
%

%#inbounds
%#realonly


[m,n]=size(C);
if m>n,
    error('Matrix cannot have more rows than columns');
end
n2=n*n;
%Put in zeros to square out maxtrix so as to not affect the score or assignment
C=[C;zeros(n-m,n)];
A=[kron(speye(n),ones(1,n));kron(ones(1,n),speye(n))];
A=A(1:end-1,:);
b=ones(2*n-1,1);
c=-C(:);
options=optimset('display','off');
[x,FVAL,EXITFLAG,OUTPUT] =linprog(c,[],[],A,b,zeros(n2,1),ones(n2,1),[],options);
X=reshape(x,n,n);
[temp,p]=max(X');
%   [c,ia,ib]=intersect(p,1:n); %Make sure p is a permutation
%   if length(ia)~=n
%      ip=setdiff(1:n,ia);
%      p(ip)=setdiff(1:n,ib);
%   end
w=-FVAL;
p=p(1:m);
x=X';
%[q,i]=unique(p);
% % Interior point solution is not a vertex!
% if length(q)~=m
%    [p,w]=munkres(-X);
% end
