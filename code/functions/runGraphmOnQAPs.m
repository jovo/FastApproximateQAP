function [mn,s1,s3,s10]=runGraphmOnQAPs(algName,problem_set)
%
%function [mn, s1,s3,s10]=runQAPs(problem_set)
% Run the QAP problem set using the graphm algorithms
%
% FIX THIS:
% -------------------------------------------------------
% FW1 is up to 30 iterations of Frank Wolf method with a single starting
%     point, the flat doubly stochastic matrix (the barycenter of the space
%     of doubly stochastic matrices.)
% FW3 is the minimum of FW1 and two starts near the barycenter of the space.
% FW10 is the minimum of FW1 and nine starts near the barycenter of the space.
% -------------------------------------------------------
%
% The input argument is one of the following:
% 'all' which will run all the QAP test problems in the QAPLIB as of
%       June 2012.
% 'path16' or 'lipa16' which will run the 16 QAPLIB test problems
%      that were run in the PATH paper and the 16 lipa[nn][a|b] problems.
% If no input or an unrecognized string is given then 'path16' 
%    is the assumed default.
% A table summarizing the results will be produced.
% Also a structure, mn containing the output of the runs including the
%     runtime, number of iteration, as well as the approximate solutions found,
%     is returned.  Also, on output  
% s1, s3, and s10 are given, which are the cost obtained using
%     1, 3, and 10 starts of the Frank Wolf method.

[a,b] = system('hostname');
if strcmp(b(1:9),'mbp.local');
    rootDir='/Users/harley/Dropbox/FastApproximateQAP';
else
    rootDir='/home/eharley/FastApproximateQAP';
end

addpath([rootDir '/code/FAQ/']) 
addpath([rootDir '/code/functions/']) 

matFileDir = [rootDir '/data/qaplib/mat_files'];

if nargin==0, 
    problem_set='';
    algName='I';
end
if iscell(problem_set)
    probs=problem_set;
else
    switch problem_set
        case'path16'
            fname=[matFileDir, '/problems16'];
            load(fname)
            probs=problems16;
        case 'lipa16'
            fname=[matFileDir, '/problems16_lip'];
            load(fname)
            probs=problems16_lip;
        case 'all'
            fname=[matFileDir, '/problems_all'];
            load(fname)
            probs=problems_all;
        otherwise
            fprintf('Default: runQAPs(''path16'')\n');
            fname=[matFileDir, '/problems16'];
            load(fname)
            probs=problems16;
    end
end

% initialize random # generator to known seed
% TODO: add this as optional argument
rand('seed',12345678);

% Loop over the starts
for i=1:10,
    display(['trial ', num2str(i)])
    % If it it the first start then start in the center of the space
    if i>1
        [mn(i).s,mn(i).p,mn(i).iters,mn(i).best,mn(i).time] = runQAPLIBProblems(probs,algName,30,-1,rootDir);
        % otherwise pick a point near the center of the space.
    else
        [mn(i).s,mn(i).p,mn(i).iters,~,mn(i).time] = runQAPLIBProblems(probs,algName,30,0,rootDir);
    end  
end

S=[mn.s]';
s1=S(1,:);
s10=min(S(1:10,:));
s3=min(S(1:3,:));
%s100=min(S)
%fprintf(1,'%12s %12s %12s %12s\n','Problem','FW1','FW2','FW10');
%for i=1:length(probs)
%    fprintf(1,'%12s %12d %12d %12d \n',probs{i},s1(i),s3(i),s10(i));
%end
% save([fname, '_results'])