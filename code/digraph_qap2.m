%% make cubic figure for QAP paper
clear all, close all
load('../data/digraph_qap.mat')


%%
clc
figure(1), clf
hold all

% figure(gcf), 
x=digraph_qap.ms;
y=digraph_qap.times;
plot(x,y,'.','markersize',8,'markeredgecolor',0.25*[1 1 1])


p=polyfit(x,y,3);


xlabel('n','fontsize',18)
ylabel('time (seconds)','fontsize',18)
title('Cubic Scaling of QAP','fontsize',18)
legend('off')
grid('on')


xs=100:100:1000;

yhat=p(1)*xs.^3+p(2)*xs.^2+p(3)*xs+p(4);
plot(xs,yhat,'k','linewidth',2)

str=['$\hat{y}=1.5\!e\!^-\!^6x^3 - 6\!e\!^-\!^4x^2 + 0.1x - 6.3$'];
text(50,900,str,'fontsize',12,'Interpreter','latex')


%
figname='../figs/digraph_qap2';
wh=[3 2]*1.5;   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-dpdf',figname)
