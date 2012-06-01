n = 279;
num_runs = 100;
num_algs = 6;

%%

%

err_chem = zeros(num_runs,num_algs);

for idx = 1:num_runs

    % Achem(the_perm,the_perm) == Achem(p_PATH,p_PATH);
    
    the_perm = dlmread(sprintf('../data/graphm/data/perm/perm_%d.txt',idx));
    the_perm = the_perm(:);
    
    iv = 1:n;
    rp = zeros(n,1);
    rp(the_perm) = iv;
    
    graphm_perms = dlmread(sprintf('../data/graphm/output/exp_out_file_%d',idx),' ',5,0);

%     p_I = graphm_perms(:,1);
%     p_rand = graphm_perms(:,5);
%     p_RANK = graphm_perms(:,3);
%     p_U = graphm_perms(:,2);
%     p_QCV = graphm_perms(:,4);
%     p_PATH = graphm_perms(:,num_algs);
    
    graphm_perms = graphm_perms(:,[1 5 3 2 4 6]);

    for ii = 1:num_algs
        err_chem(idx,ii) = sum( rp == graphm_perms(:,ii) ) / n;

        p = graphm_perms(:,ii);
        
        err_chem(idx,ii) = nnz( rp == p ) / n;
        
    end

end


%

err_gap = zeros(num_runs,num_algs);

for idx = 1:num_runs

    % Agap(the_perm,the_perm) == Agap(p_PATH,p_PATH);
    
    the_perm = dlmread(sprintf('../data/graphm/data/perm/perm_%d.txt',idx));
    the_perm = the_perm(:);
    
    iv = 1:n;
    rp = zeros(n,1);
    rp(the_perm) = iv;
    
    graphm_perms = dlmread(sprintf('../data/graphm/output/exp_out_file_gap_%d',idx),' ',5,0);

%     p_I = graphm_perms(:,1);
%     p_rand = graphm_perms(:,5);
%     p_RANK = graphm_perms(:,3);
%     p_U = graphm_perms(:,2);
%     p_QCV = graphm_perms(:,4);
%     p_PATH = graphm_perms(:,6);

    graphm_perms = graphm_perms(:,[1 5 3 2 4 6]);
    
    for ii = 1:num_algs
        
        p = graphm_perms(:,ii);
        
        err_gap(idx,ii) = nnz( rp == p ) / n;
        
    end

end

%%
figure(1);
subplot(2,1,1);
boxplot(1-err_chem,'labels',{ 'I','rand','RANK','U','QCV','PATH' },'color','k','symbol','k+');
title('Chemical')
ylabel('error (%)')
set(gca,'YTick',[0:.25:1])
ylim([0 1])

subplot(2,1,2);
boxplot(1-err_gap,'labels',{ 'I','rand','RANK','U','QCV','PATH' },'color','k','symbol','k+');
ylim([0 1])
set(gca,'YTick',[0:.25:1])
title('Electrical')
ylabel('error (%)')

wh=[3 2]*2;
figname='../figs/connectomes';
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-dpdf',figname)
print('-dpng',figname)
saveas(gcf,figname)




