%% plot qaplib results
% to generate the path16 figure, simply set bench='path16'
% to generate the lipa16 figure, simply set bench='lipa16'


clear
rootDir='~/Research/projects/primary/FastApproximateQAP/data/';

bench='path16';
switch bench
    case 'all'
        dataFileName=[rootDir, 'results/problems_all'];
    case 'path16'
        dataFileName=[rootDir, 'results/problems16'];
    case 'lipa16'
        dataFileName=[rootDir, 'results/problems16_lip'];
end

load([dataFileName, '_results'])
[foo n_qaps]=size(S);
savestuff=0;


%%

% get best known results
switch bench
    case 'all'
        best_known=mn(end).best;
        best_known(best_known<0)=NaN;
    case 'path16'
        best_known=[11156, 9896, 9504, 2298, 6194, 292, 235528, 354210, 725522, ...
            135028, 388214, 491812, 703482, 1818146, 2422002, 3139370]';
        
        path=  [18048, 19086, 16206, 5560, 8500, 300, 256320, 391270, 778284, ...
            152534, 419224, 530978, 753712, 1903872, 2555110, 3281830];
        
        qbp=   [20306, 26132, 29862, 6674, 9942, 296, 278834, 381016, 804676, ...
            165364, 455778, 550862, 799790, 1996442, 2720986, 3529402];
        % \# & Problem  &   Min   & \FAQ & \texttt{PATH} & \texttt{QBP} \\
        % 1&    chr12c &   11156 &    \textbf{13072} &   18048 	& 20306\\
        % 2&    chr15a &    9896 &    27584 &   \textbf{19086} 	& 26132\\
        % 3&    chr15c &    9504 &    17324 &   \textbf{16206} 	& 29862\\
        % 4&   chr20b &    2298 &     \textbf{3068} &    5560 		& 6674\\
        % 5&    chr22b &    6194 &    \textbf{8482} &    8500 		& 9942\\
        % 6&    esc16b &     292 &    320 &     300 		& \textbf{296}\\
        % 7&     rou12 &  235528 &    \textbf{253684} &  256320 	& 278834\\
        % 8&     rou15 &  354210 &    \textbf{371458} &  391270 	& 381016\\
        % 9&     rou20 &  725522 &    \textbf{743884} &  778284 	& 804676\\
        % 10&    tai10a &  135028 &   157954 &  \textbf{152534} 	& 165364\\
        % 11&    tai15a &  388214 &   \textbf{397376} &  419224 	& 455778\\
        % 12&    tai17a &  491812 &   \textbf{529134} &  530978 	& 550862\\
        % 13&    tai20a &  703482 &   \textbf{734276} &  753712 	& 799790\\
        % 14&    tai30a & 1818146 &  	\textbf{1894640} & 1903872 	& 1996442\\
        % 15&    tai35a & 2422002 & 	\textbf{2460940} & 2555110 	& 2720986\\
        % 16&    tai40a & 3139370 &  	\textbf{3227612} & 3281830 	& 3529402\\
        
    case 'lipa16'
        
        best_known=[
            3683
            27076
            13178
            151426
            31538
            476581
            62093
            1210244
            107218
            2520135
            169755
            4603200
            253195
            7763962
            360630
            12490441];
        
        GA=[3909, 27076, 13668, 151426, 32590, 476581, 63730, 1210244, 109809, ...
            2520135, 173172, 4603200, 258218, 7763962, 366743, 12490441];
        
        EPATH= [3885, 32081, 13577, 151426, 32247, 476581, 63339, 1210244, 109168, ...
            2520135, 172200, 4603200, 256601, 7763962, 365233, 12490441];
        %         Problem          FW1          FW2         FW10
        %         lipa20a         3791         3779         3779
        %         lipa20b        27076        27076        27076
        %         lipa30a        13571        13474        13449
        %         lipa30b       151426       151426       151426
        %         lipa40a        32109        32109        32094
        %         lipa40b       476581       476581       476581
        %         lipa50a        62962        62962        62906
        %         lipa50b      1210244      1210244      1210244
        %         lipa60a       108488       108488       108488
        %         lipa60b      2520135      2520135      2520135
        %         lipa70a       171820       171785       171611
        %         lipa70b      4603200      4603200      4603200
        %         lipa80a       256073       255779       255779
        %         lipa80b      7763962      7763962      7763962
        %         lipa90a       363937       363937       363884
        %         lipa90b     12490441     12490441     12490441
        
end

% get n_vertices
n_vertices=zeros(1,n_qaps);
for i=1:n_qaps
    n_vertices(i)=str2double(probs{i}(4:5));
    if isnan(n_vertices(i))
        n_vertices(i)=str2double(probs{i}(5:6));
    end
end
[B, IX]=sort(n_vertices);

% get times
times=zeros(10,n_qaps);
for i=1:10
    times(i,:)=mn(i).time;
end

ps=[0.05, 0.25, 0.5, 0.75, 0.95];
time_dist = quantile(times,ps);

% objective values
obj_dist = quantile(S,ps);

% error rate
err10=best_known./min(S)';


% median time sort
[B, IXtime]=sort(time_dist(3,:));

%%

figure(2),clf
switch bench
    case 'path16'
        
        % start from barycenter figure
        all=[best_known'; s1; s3; s10; s100; min([path; qbp])]';
        psoa=min([path; qbp]);
        [foo nalgs]=size(all);
        ncols=2;
        nrows=1;
        colormap('gray')
        
        subplot(ncols,nrows,1)
        alg_list=[2,nalgs];
        len=length(alg_list);
        bests=repmat(best_known,1,len);
        bar((all(:,alg_list)-bests)./repmat(best_known,1,len))
        %         set(gca,'yscale','log')
        legend('FAQ','PSOA')
        set(gca,'YTick',[0:0.5:3])
        axis('tight')
        ylabel('optimal error ratio')
        title('Performance Comparison on QAPLIB Undirected Benchmarks')
        
        
        subplot(ncols,nrows,2)
        bar(([s1]./repmat(psoa,1,1))')
        set(gca,'yscale','log')
        set(gca,'YTick',[.8, 1, 1.2])
        ylabel('FAQ / PSOA')
        xlabel('problem number')
        axis('tight')
        
        if savestuff
            figName=[rootDir, '../figs/path16'];
            wh=[3 2]*1.5;
            set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
            print('-dpdf',figName)
            print('-dpng',figName)
            saveas(gcf,figName)
        end
        
        % multiple restart figure
        figure(3), clf
        ncols=2;
        nrows=1;
        colormap('gray')
        
        all=[best_known'; s1; s3; s10; s100; min([s1; psoa])]';
        subplot(ncols,nrows,1)
        alg_list=[5,3,nalgs];
        len=length(alg_list);
        bests=repmat(best_known,1,len);
        bar((all(:,alg_list)-bests)./repmat(best_known,1,len))
        %         set(gca,'yscale','log')
        legend('FAQ_1_0_0','FAQ_3','min(FAQ,PSOA)')
        set(gca,'YTick',[0:0.25:3])
        axis('tight')
        ylabel('optimal error ratio')
        title('Performance Comparison on QAPLIB Benchmarks with Multiple Restarts')
        
        
        subplot(ncols,nrows,2)
        bar(([s100; s3; min([s1; psoa])]./repmat(min([s1; psoa]),3,1))')
        set(gca,'yscale','log')
        set(gca,'YTick',[.5:.25:2])
        ylabel('FAQ / PSOA')
        xlabel('problem number')
        axis('tight')
        
        if savestuff
            figName=[rootDir, '../figs/path16_restarts'];
            wh=[3 2]*1.5;
            set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
            print('-dpdf',figName)
            print('-dpng',figName)
            saveas(gcf,figName)
        end
        
    case 'lipa16'
        all=[best_known'; s1; s3; s10; min([GA; EPATH])]';
        psoa=min([GA; EPATH]);
        [foo nalgs]=size(all);
        ncols=2;
        nrows=1;
        colormap('gray')
        
        exclude=find(std(all,[],2)<eps);
        include=setdiff(1:16,exclude);
        
        subplot(ncols,nrows,1)
        alg_list=[2,nalgs];
        len=length(alg_list);
        bests=repmat(best_known(include),1,len);
        bar((all(include,alg_list)-bests)./bests)
        %         set(gca,'yscale','log')
        legend('FAQ','PSOA')
        %         set(gca,'YTick',[1, 1.5, 2, 2.5])
        axis('tight')
        ylabel('optimal error ratio')
        title('Performance Comparison on QAPLIB Directed Benchmarks')
        
        subplot(ncols,nrows,2)
        bar((s1(include)./repmat(psoa(include),1,1))')
        set(gca,'yscale','log')
        set(gca,'YTick',[.98:.01:1])
        ylabel('FAQ / PSOA')
        xlabel('problem #')
        axis('tight')
        
        if savestuff
            figName=[rootDir, '../figs/lipa16'];
            wh=[3 2]*1.5;
            set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
            print('-dpdf',figName)
            print('-dpng',figName)
            saveas(gcf,figName)
        end
        
        % make table
        clc
        fprintf(1,'\\hline \n');
        fprintf(1,'%12s & %12s & %12s & %12s & %12s & %12s\n','\#','Problem','Min','\texttt{FAQ}','\texttt{EPATH}','\texttt{GRAD} \\');
        fprintf(1,'\\hline \n');
        for i=1:length(probs)
            fprintf(1,'%12d & %12s & %12d & \\textbf{%12d} & %12d & %12d \\\\ \n',i,probs{i},best_known(i),s1(i),EPATH(i),GA(i));
        end
        fprintf(1,'\\hline \n');
        
        
        
    case 'all'
        % err totals
        err_tots=repmat(best_known',10,1)./S;
        
        ncols=2;
        subplot(ncols,nrows,1)
        boxplot(err_tots(:,IX),'plotstyle','compact')
        set(gca,'yscale','log')
        
        
        subplot(ncols,nrows,2)
        boxplot(err_tots(:,IXtime),'plotstyle','compact')
        set(gca,'yscale','log')
        
end