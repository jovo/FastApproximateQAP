% This tests how well we can match graphs with vertex labellings

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

% Create labellings
C = rand(n);
for i = 1:n
    % Now make it so there is a perfect labelling
    C(i,i) = 0;
end
C = C * P;

%%%%%%%%%%%%%%%%%% Test 1: alpha = 0 %%%%%%%%%%%%%%%%%%%%%%
% When alpha is set to 0 we don't care about vertex labellings so should
% find the isomorphic distance

[f, ~, sfw_P, Q] = graphm_sfw(A, B, 30, -1, C, 0);

fprintf('\nGraph matching error with alpha = 0 (optimal is 0): %g\n', f)

%%%%%%%%%%%%%%%%%% Test 2: alpha = 1 %%%%%%%%%%%%%%%%%%%%%%
% When alpha is set to 1 we only care about vertex labellings so our best
% answer should be the minimum of the LAP tr(C^T P)

[f, ~, sfw_P, Q] = graphm_sfw(A, B, 30, -1, C, 1);

% Now find the pure LAP
[p,w,x] = assign(-C);

fprintf('\nGraph matching error with alpha = 1 (optimal is %g): %g\n',-w, f)

%%%%%%%%%%%%%%%%%% Test 3: alpha = 0.5 %%%%%%%%%%%%%%%%%%%%%%
% When alpha is set to 0.5 we get a mixture between structural and labelled

[f, ~, sfw_P, Q] = graphm_sfw(A, B, 30, -1, C, 0.5);

% Now find the pure LAP
[p,w,x] = assign(-C);

fprintf('\nGraph matching error with alpha = 0.5 (optimal is 0): %g\n', f)