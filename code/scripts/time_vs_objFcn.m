clc
fnames=fieldnames(a(1).times);

markers=['+','.','o','*','x','s'];

% figure(1), clf, hold all
for j=1:10 %length(a)
    figure(j), clf, hold all
    for i=1:length(fnames)
        plot(a(j).times.(fnames{i}),-a(j).objFcns.(fnames{i}),...
            '.','marker',markers(i),...
            'markersize',15)
    end
    legend(fnames)
    xlabel('time (sec)')
    ylabel('objective function')
end

%%
absTimeDiff=zeros(length(a),length(fnames));
fracTimeDiff=absTimeDiff;
absObjFcnsDiff=absTimeDiff;
fracObjFcnsDiff=absTimeDiff;
for j=1:length(a)
    for i=1:length(fnames)
        absTime(j,i)=a(j).times.(fnames{i});
        absTimeDiff(j,i)=a(j).times.(fnames{i})-a(j).times.I;
        fracTimeDiff(j,i)=a(j).times.(fnames{i})/a(j).times.I;
        
        absObjFcn(j,i)=a(j).objFcns.(fnames{i});
        absObjFcnDiff(j,i)=a(j).objFcns.(fnames{i})-a(j).objFcns.PATH;
        fracObjFcnDiff(j,i)=a(j).objFcns.(fnames{i})/a(j).objFcns.PATH;
        
    end
end

%%

temp=10;
clc, figure(11), clf, hold all
for i=1:length(fnames)
    plot(fracTimeDiff(1:temp,i),fracObjFcnDiff(1:temp,i),...
            '.','marker',markers(i),...
            'markersize',15)
end
set(gca,'XScale','log','YScale','log')
legend(fnames,'location','best')
xlabel('time (sec)')
ylabel('objective function')

