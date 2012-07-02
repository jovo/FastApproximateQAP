clear all, clc
rootDir='~/Research/projects/primary/FastApproximateQAP/';
dataDir=[rootDir, 'data/results/'];
figDir=[rootDir, 'figs/'];
fileName= 'ErdosRenyi_results';
figName=[figDir, fileName];
load([dataDir, fileName, '.mat'])
savestuff=1;

%%
clc
figure(1), clf
hold all

% figure(gcf), 
plot(ns,times,'.','markersize',8,'markeredgecolor',0.25*[1 1 1])
xlabel('number of vertices','fontsize',18)
ylabel('time (seconds)','fontsize',18)
title('Cubic Scaling of QAP','fontsize',18)
legend('off')
grid('on')


xs=100:100:1000;

yhat=p(1)*xs.^3+p(2)*xs.^2+p(3)*xs+p(4);
plot(xs,yhat,'k','linewidth',2)

str=['$\hat{y}=6.1\times 10\!^-\!^{7}x^3 - 4.5\times 10\!^-\!^4x^2 + 0.18x - 16.7$'];
text(120,350,str,'fontsize',12,'Interpreter','latex')
axis([min(ns) max(ns) min(times) max(times)])

if savestuff
    wh=[3 2]*1.75;   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    print('-dpdf',figName)
    print('-dpng',figName)
    saveas(gcf,figName)
end