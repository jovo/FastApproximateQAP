function [Lhat ind P yhat] = graph_unlabeled_classify(As,G)
% this script classifies using a number of different approaches
% INPUT:
% As:   n x n x s array, where |V|=n, and s is the number of samples
% G:    a structure containing various fields including
%   ys: s-dimensional binary vector with class identity for each graph
%   y0: s0 dimensional vector indexing each matrix in As that is in class 0
%   y1: s1 dimensional vector indexing each matrix in As that is in class 1
%   n:  # of vertices
%   s:  # of samples
%   (optional additional fields of G)
%   case:   specifies whether in-sample, leave-one-out, or hold-out (default is in-sample)
%   nb:     if exists, do naive bayes classification
%   inc:    if exists, do incoherent classification
%   max:    if exists, do max degree signal classification
%   tru:    if exists, classify using tru signal dependent edges
%   block:  if exists, classifyin using block structure (not working)
%   Nmax:   # vertices to select when doing max degree
%   Ninc:   # edges to select when doing incoherent classifier
% 
% OUTPUT:
%   Lhat:   misclassification rate for each algorithm implemented
%   ind:    index of edges used for each algorithm implemented
%   P:      parameter estimates
%   yhat:   list of estimated class identity for each matrix

if ~isfield(G,'case'), G.case='in'; end

if strcmp(G.case,'in'), P = get_params(As,G); H=G; end % in-sample classifier for debugging purposes


for i=1:G.s

    Bs=As;
    
    % set up stuff for leave-one-out classification
    if strcmp(G.case,'loo'), 
        H = G;
        H.ys(i) = [];
        H.y0(H.y0==i) = [];
        H.y1(H.y1==i) = [];
        
        A=As(:,:,i);
        
        for j=[H.y0 H.y1]
            B=Bs(:,:,j);
            [f,fwq,x,iter]=sfw(A,-B);
            Bs(:,:,j)=B(fwq,fwq);
        end
        
        P = get_params(Bs,H);     % update parameters using only other samples
    end

    if isfield(G,'tru')             % classify using only true signal edges
        yhat.tru(i)  = ie_classify(Bs(:,:,i),P,G.tru);
    end

    if isfield(G,'nb')              % naive bayes classifier
        ind(i).nb   = 1:G.n^2;      % list of edges to use
        yhat.nb(i)  = ie_classify(Bs(:,:,i),P,ind(i).nb); % estimated class identities
    end

    if isfield(G,'inc')             % incoherent edge classifier
        ind(i).inc = get_inc_edges(P.d_opt,H);
        yhat.inc(i)= ie_classify(Bs(:,:,i),P,ind(i).inc);
    end

    if isfield(G,'max')             % max degree signal independent edge classifier
        ind(i).max  = get_max_edges(P.d_opt,H);
        yhat.max(i) = ie_classify(Bs(:,:,i),P,ind(i).max);
    end

    if isfield(G,'block')           % stochastic blockmodel classifier
        B = get_blocks(As(:,:,i),G);
        yhat.block(i) = MR_classify_var1(B,P,1:G.b);
    end

end

fn=fieldnames(yhat);                % names of classifiers
for i=1:length(fn)                  % for each, compute Lhat
    Lhat.(fn{i})=sum(abs(yhat.(fn{i})-G.ys))/G.s;
end
if isfield(G,'block'), Lhat.block = sum(abs(yhat.block-G.ys))/G.s; end

end

function y = ie_classify(datum,P,ind)
% this function classifies using independent edge assumption.  each edge
% could be distribution according to a poisson distribution, or a
% bernoulli. the class conditional posterior is computed as appropriate.
% 
% INPUT
% datum:    the graph to be classified
% P:        structure of parameter estimates to use for classification
% ind:      indices to use in classifier
% 
% OUTPUT:
% y:        estimated class

if any(datum(ind))>1            % if poisson
    post0=sum(sum(datum(ind).*P.lnE0(ind) - P.E0(ind)));
    post1=sum(sum(datum(ind).*P.lnE1(ind) - P.E1(ind)));
else                            % if bernoulli
    post0=sum(datum(ind).*P.lnE0(ind)+(1-datum(ind)).*P.ln1E0(ind));
    post1=sum(datum(ind).*P.lnE1(ind)+(1-datum(ind)).*P.ln1E1(ind));
end

[foo bar] = sort([post0, post1]);
y=bar(2)-1;


end



