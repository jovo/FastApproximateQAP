clear, clc

n   = 10;
s   = 1000;
p0  = 0.5;
p1  = 0.5;


% make egg
m       = 3;            % size of egg
thesem  = 1:m;          % which vertices are the signal containing vertices
p0_egg  = 0.25;
p1_egg  = 0.75;

% make E0
E0          = p0*ones(n); 
E0(thesem,thesem) = p0_egg;

% make E1
E1          = p1*ones(n); 
E1(thesem,thesem) = p1_egg;

% generate data
As  = zeros(n,n,s);
As(:,:,1:s/2)   = repmat(E0,[1 1 s/2]) > rand(n,n,s/2);
As(:,:,s/2+1:s) = repmat(E1,[1 1 s/2]) > rand(n,n,s/2);

% make data unlabeled
for i=1:s
    q=randperm(n);
    A=As(:,:,i);
    As(:,:,i)=A(q,q);
end

ys = [zeros(1,s/2) ones(1,s/2)];
G = get_constants(As,ys);

G.tru=[]; G.nb=[]; G.inc=[]; G.max=[]; G.Ninc=m^2; G.Nmax=m;

[Lhat ind P yhat] = graph_unlabeled_classify(As,G)

