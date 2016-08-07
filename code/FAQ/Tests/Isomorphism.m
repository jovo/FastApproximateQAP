% FAQ should solve QAP exactly when A and B are adjacency matrices  of
% simple graphs and are isomorphic to one another

% Define test dimension
p = 5;

% Define A as a random adjacenecy matrix of a simple graph
% To do this we need p(p-1)/2 Bernoulli trials for each edge
A = zeros(p);
for i = 1:p-1
    for j = i+1:p
        A(i,j) = rand() > 0.5;
        A(j,i) = A(i,j);
    end
end

% Now generate a random permutation matrix
P = perm2mat( randperm(p) );

% Produce an isomorphic B
B = P * A * P';

% Now run sfw and check the permutation transforms A onto B
[~, sfw_p] = sfw(-A, B);  % Note we use -A as sfw solves QAP
sfw_P = perm2mat(sfw_p);

assert(all(all(A == sfw_P * B * sfw_P')))

