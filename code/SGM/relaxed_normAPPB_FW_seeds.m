% This is an implementation of the Frank-Wolfe method for solving the following optimization problem:
% min_{P \in D} ||AP - PB||_F^2
% where D is the set of doubly stochastic matrices, A and B are square matrices.
% This optimization problem comes from the Graph Matching Problem, when relaxing the set of permutation matrices to its convex hull, which is the set of doubly stochastic matrices. Here A and B are de adjacency matrices of the two graphs to match. 
% When, in addition, the first nodes of each graph are already in correspondance (called "seeds"), this information can be added to help the algorithm. 
% Here, the argument "seeds" is an integer which indicates the number of nodes that are already in correspondance. It's assumed that these nodes are the first ones.
%
%
% This code was used in the following publications:
%
% On spectral properties for graph matching and graph isomorphism problems. Marcelo Fiori and Guillermo Sapiro. Information and Inference (2015) 4 (1): 63-76
%
% Graph Matching: Relax at Your Own Risk. Vince Lyzinski, Donniell Fishkind, Marcelo Fiori, Joshua T. Vogelstein, Carey E. Priebe, and Guillermo Sapiro. IEEE Transactions on Pattern Analysis and Machine Intelligence Vol. 38, 1, pp. 60-73, 2016
%
% Robust Multimodal Graph Matching: Sparse Coding Meets Graph Matching. Marcelo Fiori, Pablo Sprechmann, Joshua Vogelstein, Pablo Musé, Guillermo Sapiro. Advances in Neural Information Processing Systems 26 (NIPS 2013). 
%
%
% Written by Marcelo Fiori, 2013.
% Redistribution and use in source and binary forms, with or without modification, are permitted provided that the redistribution must retain the above copyright notice, and this list of conditions.


function [P,Pp]=relaxed_normAPPB_FW_seeds(A,B,seeds)

AtA = A'*A;
BBt = B*B';
p=size(A,1);

f1 = @(P) norm(A*P-P*B,'fro')^2;

tol=5e-2;
tol2=1e-5;

P=ones(p)/(p-seeds);
P(1:seeds,1:seeds)=eye(seeds);

f=f1(P);
var=1;

while ~(abs(f)<tol) && (var > tol2)
    fold=f;

    grad = (AtA*P -A'*P*B - A*P*B' + P*BBt);
    
    grad(1:seeds,:)=0;
    grad(:,1:seeds)=0;
    G=round(grad);
    corr=lapjv(G(seeds+1:end,seeds+1:end),0.01);  
    Ps=perm2mat(corr);
    
    Ps1=eye(p);
    Ps1(seeds+1:end,seeds+1:end) = Ps;
    Ps=Ps1;
    
    C = A*(P-Ps) + (Ps-P)*B; 
    D = A*Ps-Ps*B;
    
    aq = trace(C*C');
    bq = trace(C*D'+D*C');
    aopt = -bq/(2*aq);

    Ps4 = aopt*P + (1-aopt)*Ps;
    
    f=f1(Ps4);
    P=Ps4;
    
    var=abs(f-fold);

end

    corr=lapjv(-P,0.01);
    Pp=perm2mat(corr);

end
