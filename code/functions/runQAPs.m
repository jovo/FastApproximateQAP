function [mn,s1,s3,s10]=runQAPs(problem_set)
%function [mn, s1,s3,s10]=runQAPs(problem_set)
% Run the QAP problem set using FW1, FW3, and FW10
% where
% FW1 is up to 30 iterations of Frank Wolf method with a single starting
%     point, the flat doubly stochastic matrix (the barycenter of the space
%     of doubly stochastic matrices.)
% FW3 is the minimum of FW1 and two starts near the barycenter of the space.
% FW10 is the minimum of FW1 and nine starts near the barycenter of the space.
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
% John M. Conroy, June 2012
% IDA Center for Computing Sciences
% conroy@super.org
%  (c) 2012, Institute for Defense Analyses, 4850 Mark Center Drive, 
%  Alexandria, Virginia, 22311-1882; 703-845-2500.
%  This material may be reproduced by or for the U.S. 
%  Government pursuant to the copyright license under 
%  the clauses at DFARS 252.227-7013 and 252.227-7014.
%
% To run this code, one must be in the /qaplib directory, with the /code/
% directory in the path.

rootDir='~/Research/projects/primary/FastApproximateQAP/qaplib/';
cd(rootDir)
path('../code/FAQ/',path) 
path('../code/other/',path) 


if nargin==0, 
    problem_set='';
end
if iscell(problem_set)
    probs=problem_set;
else
    switch problem_set
        case'path16'
            fname=[rootDir, 'mat_files/problems16'];
            load(fname)
            probs=problems16;
        case 'lipa16'
            fname=[rootDir, 'mat_files/problems16_lip'];
            load(fname)
            probs=problems16_lip;
        case 'all'
            fname=[rootDir, 'mat_files/problems_all'];
            load(fname)
            probs=problems_all;
        otherwise
            fprintf('Default: runQAPs(''path16'')\n');
            fname=[rootDir, 'mat_files/problems16'];
            load(fname)
            probs=problems16;
    end
end
rand('seed',12345678);
% Loop over the starts
for i=1:10,
    display(['trial ', num2str(i)])
    % If it it the first start then start in the center of the space
    if i>1
        [mn(i).s,mn(i).p,mn(i).iters,mn(i).best,mn(i).time]=run_qaplib_problems(probs,30,-1);
        % otherwise pick a point near the center of the space.
    else
        [mn(i).s,mn(i).p,mn(i).iters,~,mn(i).time]=run_qaplib_problems(probs,30);
    end
    
end;
S=[mn.s]';
s1=S(1,:);
s10=min(S(1:10,:));
s3=min(S(1:3,:));
%s100=min(S)
fprintf(1,'%12s %12s %12s %12s\n','Problem','FW1','FW2','FW10');
for i=1:length(probs)
    fprintf(1,'%12s %12d %12d %12d \n',probs{i},s1(i),s3(i),s10(i));
end
% save([fname, '_results'])
