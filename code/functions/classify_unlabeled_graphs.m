function [LAP QAP] = classify_unlabeled_graphs(Atrn,Atst,ytst,P,alg)
% this function classifies unlabeled graphs, in the following steps:
%
% (i) permute Atst to be like both a single Atrn in class 0 and a single
% Atrn in class 1
% (ii) compute the likelihood of the permuted Atsts coming from either of
% the classes
%
% we approximate solving the permuation by using LAP and/or QAP
% initial conditions can be set for both
% QAP is iterative, so we can set the number of iterations
%
% INPUT:
%   Atrn:   an n x n x S/2 array of adjacency matrices for training
%   Atst:   an n x n x S/2 array of adjacency matrices for testing
%   ytst:   list of classes for test data
%   P:      structure of parameters necessary for naive bayes classification
%   alg:    a structure of settings for the algorithms
%
% OUTPUT:
%   LAP:    a structure containing all the LAP results
%   QAP:    a structure containing all the QAP results

%% initialize stuff

siz = size(Atrn);
n   = siz(1);       % # of vertices
S   = 2*siz(end);   % total # of samples

S1  = sum(ytst);    % # of test samples from class 1
S0  = S-S1;         % # of test samples from clsss 0

if alg.truth_start == false     % if not starting at the truth, generate permuations for the test graphs
    tst_ind=zeros(S/2,n);
    for j=1:S/2
        tst_ind(j,:) = randperm(n);
    end
end

if strcmp(alg.names(1),'LAP'),
    LAP.time    = 0;            
    LAP.do      = true;
    LAP.yhat    = NaN(S/2,1);
    LAP.correct = NaN(S/2,1);
    LAP.work0   = NaN(S/2,1);
    LAP.work1   = NaN(S/2,1);
else LAP.do = false;
end

if strcmp(alg.names(2),'QAP'),
    QAP.time    = 0;
    QAP.do      = true;
    QAP.yhat    = NaN(S/2,alg.QAP_max_iters);
    QAP.correct = NaN(S/2,alg.QAP_max_iters);
    QAP.work0   = NaN(S/2,alg.QAP_max_iters);
    QAP.work1   = NaN(S/2,alg.QAP_max_iters);
    QAP.obj0    = NaN(S/2,alg.QAP_max_iters);
    QAP.obj1    = NaN(S/2,alg.QAP_max_iters);
else QAP.do = false;
end

for j=1:S/2
    
    if alg.truth_start == false
        A = Atst(tst_ind(j,:),tst_ind(j,:),j);
    else
        A = Atst(:,:,j);
    end
    
    % do LAP
    if LAP.do
        starttime = cputime;
        LAP_ind0 = munkres(-Atrn(:,:,1)*A');     %Note that sum(sum(-A(LAP,:).*Atst))==f
        LAP.work0(j)=norm(Atrn(:,:,1)-A(LAP_ind0,:)) < norm(Atrn(:,:,1)-A); % check that LAP work
        
        LAP_ind1 = munkres(-Atrn(:,:,2)*A');     %Note that sum(sum(-A(LAP,:).*A))==f
        LAP.work1(j)=norm(Atrn(:,:,2)-A(LAP_ind1,:)) < norm(Atrn(:,:,2)-A); % check that LAP work
        
        LAP_lik0=sum(sum(A(LAP_ind0,:).*P.lnE0+(1-A(LAP_ind0,:)).*P.ln1E0));
        LAP_lik1=sum(sum(A(LAP_ind1,:).*P.lnE1+(1-A(LAP_ind1,:)).*P.ln1E1));
        
        [~, bar] = sort([LAP_lik0, LAP_lik1]);
        LAP.yhat(j)=bar(2)-1;
        LAP.correct(j)=(LAP.yhat(j)==ytst(j));
        LAP.time = LAP.time + cputime-starttime;
    end
    
    % do QAP
    if QAP.do
        
        starttime = cputime;
        [~,~,~,iter0,~,QAP_inds0]=sfw(Atrn(:,:,1),-A,alg.QAP_max_iters,alg.QAP_init);
        
        QAP.obj0(j,1) = norm(Atrn(:,:,1)-A);
        for ii=1:iter0
            QAP.obj0(j,ii+1) = norm(Atrn(:,:,1) - A(QAP_inds0{ii},QAP_inds0{ii}));
            QAP.work0(j,ii)= QAP.obj0(j,ii+1) < QAP.obj0(j,1); % check that QAP works
        end
        
        [~,~,~,iter1,~,QAP_inds1]=sfw(Atrn(:,:,2),-A,alg.QAP_max_iters,alg.QAP_init);
        QAP.obj1(j,1)  = norm(Atrn(:,:,1)-A);
        for ii=1:iter1
            QAP.obj1(j,ii+1) = norm(Atrn(:,:,2) - A(QAP_inds1{ii},QAP_inds1{ii}));
            QAP.work1(j,ii)= QAP.obj1(j,ii+1) < QAP.obj1(j,1); % check that QAP works
        end
        
        for ii=1:min(iter0,iter1)
            QAP_lik0=sum(sum(A(QAP_inds0{ii},QAP_inds0{ii}).*P.lnE0+(1-A(QAP_inds0{ii},QAP_inds0{ii})).*P.ln1E0));
            QAP_lik1=sum(sum(A(QAP_inds1{ii},QAP_inds1{ii}).*P.lnE1+(1-A(QAP_inds1{ii},QAP_inds1{ii})).*P.ln1E1));
            
            [~, bar] = sort([QAP_lik0, QAP_lik1]);
            QAP.yhat(j,ii)=bar(2)-1;
            QAP.correct(j,ii)=(QAP.yhat(j)==ytst(j));
        end
        QAP.time = QAP.time + cputime-starttime;
        
    end
end

if LAP.do
    LAP.Lhat = 1-mean(LAP.correct);
    LAP.Lvar = var(LAP.correct);
end

if QAP.do
    for ii=1:alg.QAP_max_iters
        corrects = QAP.correct(:,ii);
        keeper   = ~isnan(corrects);
        corrects = corrects(keeper);
        
        QAP.Lhat(ii)= 1-mean(corrects);
        QAP.Lvar(ii)= var(corrects);
        QAP.Lstd(ii)= std(corrects);
        QAP.num(ii) = length(corrects);
        QAP.obj0_avg(ii+1) = mean(QAP.obj0(keeper,ii+1));
        QAP.obj1_avg(ii+1) = mean(QAP.obj1(keeper,ii+1));
        QAP.obj0_var(ii+1) = var(QAP.obj0(keeper,ii+1));
        QAP.obj1_var(ii+1) = var(QAP.obj1(keeper,ii+1));
    end
    QAP.obj0_avg(+1) = mean(QAP.obj0(:,+1));
    QAP.obj1_avg(+1) = mean(QAP.obj1(:,+1));
    QAP.obj0_var(+1) = var(QAP.obj0(:,+1));
    QAP.obj1_var(+1) = var(QAP.obj1(:,+1));
end
