results=[    11156    12176    13072    13072    13072    18048
    9896     9896    17272    17272    27584    19086
    9504    10960    14274    14274    17324    16206
    2298     2786     3068     3068     3068     5560
    6194     7218     7876     7876     8482     8500
    292      292      294      294      320      296
    235528   235528   238134   253684   253684   256320
    354210   356654   371458   371458   371458   381016
    725522   730614   743884   743884   743884   778284
    135028   135828   148970   157954   157954   152534
    388214   391522   397376   397376   397376   419224
    491812   496598   511574   511574   529134   530978
    703482   711840   721540   721540   734276   753712
    1818146  1844636  1890738  1894640  1894640  1903872
    2422002  2454292  2460940  2460940  2460940  2555110
    3139370  3187738  3194826  3194826  3227612  3281830];


rel_results=results./repmat(results(:,1),1,6);

err=rel_results(:,[2 5 6])-1;

%%


%%
clc
figure(1), clf,
fs=12;
bar(err)
colormap('gray')
set(gca,'YScale','log')
hold on
plot(ix-.2,[1 1 1],'.','MarkerSize',16,'MarkerFaceColor','k','MarkerEdgeColor','k')

ylabel('Error','FontSize',fs)
title('Approximate QAP Performance on QAP Benchmark Library','fontsize',fs)
legend('QAP_1_0_0','QAP_1','PSOA','Location','best')
legend1 = legend(gca,'show');
set(legend1,'Orientation','horizontal','Location','NorthEast','LineWidth',0,...
    'FontName','Courier','FontSize',8);

set(gca,'FontSize',fs)
% axis([0 17 10 -100])

set(gca,'XTick',1:16,'XTickLabel',[...
    {'chr12c'}; ...
    {'chr15a'}; ...
    {'chr15c'}; ...
    {'chr20b'}; ...
    {'chr22b'}; ...
    {'esc16b'}; ...
    {'rou12'}; ...
    {'rou15'}; ...
    {'rou20'}; ...
    {'tai10a'}; ...
    {'tai15a'}; ...
    {'tai17a'}; ...
    {'tai20a'}; ...
    {'tai30a'}; ...
    {'tai35a'}; ...
    {'tai40a'};],...
    'FontSize',fs);

xticklabel_rotate;
set(gca,'YLim',[0.004 3])

wh=[6 3]*1.2;
figname='../../figs/benchmarks';
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-dpdf',figname)
print('-dpng',figname)



