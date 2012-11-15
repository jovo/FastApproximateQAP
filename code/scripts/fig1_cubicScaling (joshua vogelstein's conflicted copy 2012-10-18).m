clear all, clc
rootDir='~/Research/projects/primary/FastApproximateQAP/';
dataDir=[rootDir, 'data/results/'];
figDir=[rootDir, 'figs/'];
fileName= 'ErdosRenyi_results';
figName=[figDir, fileName];
load([dataDir, fileName, '.mat'])
savestuff=0;

%%
clc
figure(1), clf
hold all

% figure(gcf), 
plot(ns,times,'.','markersize',20,'markeredgecolor',0*[1 1 1])
xlabel('number of vertices','fontsize',18)
ylabel('time (seconds)','fontsize',18)
title('Cubic Scaling of QAP','fontsize',18)
legend('off')
grid('on')

xs=100:100:1000;

yhat=p(1)*xs.^3+p(2)*xs.^2+p(3)*xs+p(4);
plot(xs,yhat,'k','linewidth',2,'color',0.75*[1 1 1])

formatSpec = '%0.2G';
str1=['$\hat{y}=3.4 \!\times \! 10^{-9} \, x^3 + 2.7 \!\times\! 10^{-5}\, x^2$'];
str2=['$- 0.0026 x + 0.25$'];
fs=14;
% , ...
%     num2str(p(1), formatSpec),'x^3 - ', ...
%     num2str(p(2), formatSpec), 'x^2 + ', ...
%     num2str(p(3), formatSpec), 'x - ', ...
%     num2str(p(4), formatSpec), '$'];
text(120,300,str1,'fontsize',fs,'Interpreter','latex')
text(200,265,str2,'fontsize',fs,'Interpreter','latex')
axis([min(ns) max(ns) min(times) max(times)])
set(gca,'fontsize',fs)

if savestuff
    wh=[3 2]*1.75;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print('-dpdf',figName)
    print('-dpng',figName)
    saveas(gcf,figName)
end