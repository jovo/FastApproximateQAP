function [Lhat Lsem yhat] = naive_bayes_classify(preds,classes,params,subspace)
% this script does naive bayes classification, using a (sub-)space
% specified by user
%
% INPUT:
%   preds:    an array (currently typically n x n x s) of predictor semiables
%   classes:    a list of classes (typically 1 x s)
%   params:   a structure containing parameters for the algorithm to use
%   subspace: structure containing name and index of each subspace used
%
% OUTPUT:
%   Lhat:   misclassification rate for each algorithm implemented
%   Lsem:   misclassification standard error of the mean for each algorithm implemented
%   yhat:   list of estimated class identity for each graph

siz = size(preds);
n   = siz(1);                           % # of semiables
S   = siz(end);                         % # of samples (robust to preds being an arbitrary array)
sqrtS = sqrt(S);
if nargin == 3, subspace(1).all=1:n^2; end

poiss       = any(preds(:))>1;          % whether data is Poisson or Bernoulli

fn = fieldnames(subspace);              % names of subspaces
n_subspace = length(fn);                % # of different subspace

for i=1:S
    datum = preds(:,:,i);               % this line only makes sense when data are graphs    
    for j=1:n_subspace
        ind = subspace.(fn{j});         % for code legibility
        
        data_tmp=datum(ind);
        
        if poiss                        % if poisson
            post0=sum(sum(data_tmp.*params.lnE0(ind) - params.E0(ind)))+params.lnprior0;
            post1=sum(sum(data_tmp.*params.lnE1(ind) - params.E1(ind)))+params.lnprior1;
        else                            % if bernoulli
            post0=sum(data_tmp.*params.lnE0(ind)+(1-data_tmp).*params.ln1E0(ind))+params.lnprior0;
            post1=sum(data_tmp.*params.lnE1(ind)+(1-data_tmp).*params.ln1E1(ind))+params.lnprior1;
        end
        
        [~, bar] = sort([post0, post1]); % find the bigger one
        yhat.(fn{j})(i)=bar(2)-1;
    end
end

siz=size(classes);
if siz(1)>siz(2), classes=classes'; end

for j=1:n_subspace                              % for each, compute stats
    correct_vect = abs(yhat.(fn{j})-classes);
    Lhat.(fn{j}) = mean(correct_vect);          % percent correct
    Lsem.(fn{j}) = std(correct_vect)/sqrtS;     % s.e.m of correct
end