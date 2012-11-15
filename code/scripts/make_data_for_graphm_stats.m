% A : vector of structures
% 
% A(i).name       name of problem
% A(i).sol_p      permutation of optimal solution, [] if doesn't exist
% A(i).obj_p      objective function value of optimal solution, nan if doesn't exist
% 
%     A(i).FAQ.p         solution
%     A(i).FAQ.t         minimum of run times (.timings)
%     A(i).FAQ.laps      LAP count
%     A(i).FAQ.obj       objective function value
% 
%     A(i).FAQ.timings    vector of run times
% 
% FAQ
% I
% PATH
% QCV
% rand
% RANK
% U



%%

addpath('/Users/harley/Dropbox/FastApproximateQAP/code/functions');

load('/Users/harley/Dropbox/FastApproximateQAP/data/qaplib/mat_files/problems_all.mat');

cd('/Users/harley/Dropbox/FastApproximateQAP/data/graphm/09.20.2012')
load FAQ_20-Sep-2012_20.01.00.mat
mn = ans;

FAQ  = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);
I    = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);
PATH = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);
QCV  = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);
rand = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);
RANK = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);
U    = struct('p',[], 't',nan, 'laps',nan, 'obj',nan, 'timings',[]);

A = struct('name','','sol_p',[],'obj_p',nan,'FAQ',FAQ, 'I',I, 'PATH',PATH, 'QCV',QCV, 'rand',rand, 'RANK',RANK, 'U',U);

A(137) = A;

for j = 1:137

    A(j).name = problems_all{j};

    [~,~,p,s]=qap_read(A(j).name,'/Users/harley/Dropbox/FastApproximateQAP/data/qaplib');
	A(j).sol_p = p;

    A(j).obj_p = s;
end


% FAQ
load FAQ_20-Sep-2012_20.01.00.mat
mn = ans;
i = 2;
for j = 1:137
	A(j).FAQ.obj = mn(i).s(j);
	A(j).FAQ.p = mn(i).p{j};
	A(j).FAQ.laps = mn(i).iters(j);
    
    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
    A(j).FAQ.t = min(x);
        
end

% I
load I_14-Sep-2012_16.00.49.mat
i = 2;
for j = 1:137
	A(j).I.obj = mn(i).s(j);
	A(j).I.p = mn(i).p{j};
	A(j).I.laps = mn(i).iters(j);

    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
    A(j).I.t = min(x);

end

% PATH
load PATH_14-Sep-2012_19.48.20.mat
i = 2;
for j = 1:137
	A(j).PATH.obj = mn(i).s(j);
	A(j).PATH.p = mn(i).p{j};
	A(j).PATH.laps = mn(i).iters(j);

    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
	A(j).PATH.t = min(x);

end

% QCV
load QCV_14-Sep-2012_16.06.42.mat
i = 2;
for j = 1:137
	A(j).QCV.obj = mn(i).s(j);
	A(j).QCV.p = mn(i).p{j};
	A(j).QCV.laps = mn(i).iters(j);

    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
	A(j).QCV.t = min(x);

end

% rand
load rand_14-Sep-2012_16.07.23.mat
i = 2;
for j = 1:137
	A(j).rand.obj = mn(i).s(j);
	A(j).rand.p = mn(i).p{j};
	A(j).rand.laps = mn(i).iters(j);

    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
    A(j).rand.t = min(x);

end

% RANK
load RANK_14-Sep-2012_16.02.25.mat
i = 2;
for j = 1:137
	A(j).RANK.obj = mn(i).s(j);
	A(j).RANK.p = mn(i).p{j};
	A(j).RANK.laps = mn(i).iters(j);
    
    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
	A(j).RANK.t = min(x);

end

% U
load U_14-Sep-2012_16.01.37.mat
i = 2;
for j = 1:137
	A(j).U.obj = mn(i).s(j);
	A(j).U.p = mn(i).p{j};
	A(j).U.laps = mn(i).iters(j);
    
    x = [];
    for k = 1:10
        x(k) = mn(k).time(j);
    end
	A(j).U.t = min(x);

end

%%

for j = 1:137
    A(j).PATH.laps = mn_PATH.iters(j);
end

%%

% Check the timing data to make sure it's not nuts
% A(1).PATH.t

% 
% for i = 1:137
%     xx = A(i).FAQ.t;
%     if xx < 0 || xx > 100
%         i
%     end
%     
%     xx = A(i).I.t;
%     if xx < 0 || xx > 100
%         i
%     end
%     
%     xx = A(i).PATH.t;
%     if xx < 0 || xx > 100
%         i
%     end
%     
%     xx = A(i).QCV.t;
%     if xx < 0 || xx > 100
%         i
%     end
%     
%     xx = A(i).rand.t;
%     if xx < 0 || xx > 100
%         i
%     end
%     
%     xx = A(i).RANK.t;
%     if xx < 0 || xx > 100
%         i
%     end
%     
%     xx = A(i).U.t;
%     if xx < 0 || xx > 100
%         i
%     end
% end
% 
% 
