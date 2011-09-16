%% make cubic figure for QAP paper
clear all, close all
load('../data/digraph_qap.mat')


%%

figure(1)
hold all

% figure(gcf), 
x=digraph_qap.ms;
y=digraph_qap.times;
plot(x,y,'.','markersize',8,'markeredgecolor',0.25*[1 1 1])


p=polyval(x,y,3);


xlabel('n','fontsize',18)
ylabel('time (seconds)','fontsize',18)
legend('off')
grid('on')

str=['y=ax+b'];
text(0,100,str)


figname='digraph_qap2';
% print('-dpdf',figname)
