% FAQ test where A and B are adjacency matrices  of
% simple graphs and are isomorphic to one another. We hope to see the graph
% matching distance as 0 although this will not be the case if rQAP
% converges to a local minimum.

% Define test dimension
n = 10;

% Define A as a random adjacenecy matrix of a simple graph
% To do this we need p(p-1)/2 Bernoulli trials for each edge
A = zeros(n);
for i = 1:n-1
    for j = i+1:n
        A(i,j) = rand() > 0.5;
        A(j,i) = A(i,j);
    end
end

% Now generate a random permutation matrix
P = perm2mat( randperm(n) );

% Produce an isomorphic B
B = P * A * P';

% Now run sfw and check the permutation transforms A onto B
% [~, sfw_p] = sfw(-A, B);  % Note we use -A as sfw solves QAP
% sfw_P = perm2mat(sfw_p);

[f, ~, sfw_P, Q] = graphm_sfw(A, B);

fprintf('\nGraph matching error (optimal is 0): %g\n', f)


