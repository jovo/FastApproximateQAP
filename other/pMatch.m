function [ X,permutation ] = pMatch( A,B )

% type the following commands in MATLAB:
%
% n=40;
% A=round(rand(n,n));A=A-triu(A);A=A+A';p1=randperm(n);B=A(p1,p1);
% [X,p2]=pMatch(A,B)
%
% This generates a random n-by-n ErdosRenyi(1/2)graph adjacency
% matrix A, and random permutation p1, and p1-permuted-version 
% of A called B.
% X is a doubly stochastic matrix solving AX=XB, and if X is integral
% then p2 is the reconstructed permutation (hopefully is p1?)

[m,n]=size(A);
M=[kron(eye(n),A)-kron(B',eye(n));kron(eye(n),ones(1,n));kron(ones(1,n),eye(n))];
b=sparse([zeros(n^2,1);ones(n,1);ones(n,1)]);
H=sparse([M eye(n^2+2*n)]);
c=sparse([zeros(n^2,1);ones(n^2+2*n,1)]);
startx=sparse([zeros(n^2,1); b]);
lb=sparse(zeros(2*n^2+2*n,1));
y=linprog(c,[],[],H,b,lb,[],startx);
X=zeros(n,n);
for i=1:n
    X(:,i)=y((i-1)*n+1:i*n);
end
permutation=inf*ones(1,n);
for i=1:n
    for j=1:n
        if abs(X(i,j)-1)<10^(-2)
            permutation(j)=i;
        end
    end
end


