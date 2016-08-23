function [f,p,P,Q,iter,fs,myps]=graphm_sfw(A,B,IMAX,x0,C,alpha)
% A wrapper around the sfw function for solving QAP to return values
% instead relevant to graph matching

% Perform at most IMAX iterations of the Frank-Wolfe method to compute an
% approximate solution to the graph matching problem given the
% matrices A and B. A and B should be square and the same size.  The method
% seeks a permutatation p with corresponding matrix P(p) which minimizes
%       g(p)=|| A - P(p)^T B P(p) ||_Fro
% by minimizing
%       f(p)=sum(sum(A.*B(p,p)))
% Convergence is declared if a fix point is encountered or if the projected
% gradient has 2-norm of 1.0e-4 or less.
% IMAX is optional with a default value of 30 iterations.
%     If IMAX is set to 0.5 then one iteration of FW is performed with no
%     line search.  This is Carey Priebe's LAP approximation to the QAP.
% The starting point is optional as well and its default value is
% ones(n)/n, the flat doubly stochastic matrix.
% x0 may also be
%   -1 which signifies a "random" starting point should be used.
%       here the start is given by
%       0.5*ones(n)/n+sink(rand(n),10)
%       where sink(rand(n),10) performs 10 iterations of Sinkhorn balancing
%       on a matrix whose entries are drawn from the uniform distribution
%       on [0,1].
%   x0 may also be a user specified n by n doubly stochastic matrix.
%   x0 may be a permutation vector of size n.
% On output:
%     f=|| A - P(p)^T B P(p) ||_Fro, where
%     p is the permutation found by FW after projecting the interior point
%         to the boundary.
%     P is the permutation matrix corresponding to permutation p
%     Q is the doubly stochastic matrix (interior point) computed by the FW
%       method
%     iter is the number of iterations of FW performed.
%     fs is the list of fs for each iteration
%     myps is the list of myps for each iteration

    if ~exist('IMAX','var')
        [~,p,Q,iter,fs,myps] = sfw(-A, B);
        alpha = 1;
        C = zeros(size(A));
    elseif ~exist('x0','var')
        [~,p,Q,iter,fs,myps] = sfw(-A, B, IMAX);
        alpha = 1;
        C = zeros(size(A));
    elseif ~exist('C','var')
        [~,p,Q,iter,fs,myps] = sfw(-A, B, IMAX, x0);
        alpha = 1;
        C = zeros(size(A));
    elseif ~exist('alpha','var')
        [~,p,Q,iter,fs,myps] = sfw(-A, B, IMAX, x0, C);
        alpha = 1;
    else
        [~,p,Q,iter,fs,myps] = sfw(-A, B, IMAX, x0, C, alpha);
    end
    
    P = perm2mat(p);
    P = P';
    f = (1-alpha)*norm(A - P' * B * P, 'fro') + alpha*trace(C' * P);

end

