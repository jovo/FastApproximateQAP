function [A,B,p,s]=qap_read(problem,directory)
%function [A,B,p,s]=qap_read(problem,directory)
% Read a qap problem and best solution
% for the specified problem named "problem"
% Files are assumed to be in qaplib format
% with the problem file in the subdirectory
% qaplib/qapprob and the solution file if
% it exists in qaplib/qapsoln, unless the
% second input argument "directory" is present.
% In this case, qaplib is replaced by the
% value of directory.
% If a solution file exists for problem then
% p will be the permutation solution vector giving
% the best know solution, which has a value of
% s=sum(sum(A.*B(p,p)));
% Otherwise, p and s will be set to the null matrix.
% John M. Conroy, IDA/CCS
% April 1999

if exist('directory','var')==0
   directory='qaplib';
end
fp=fopen([directory,'/qapprob/',problem,'.dat'],'r');
if fp==-1
   disp(['Cannot open file:',directory,'/qapprob/',problem,'.dat']); 
   return
end
n=fscanf(fp,'%d',1);
A=fscanf(fp,'%f',[n,n])';
B=fscanf(fp,'%f',[n,n])';
fclose(fp);
% If a file with the best known solution is present, read it in.
fp=fopen([directory,'/qapsoln/',problem,'.sln'],'r');
if fp==-1
   p=[];
   s=-1;
   return
end
n=fscanf(fp,'%d',1);
s=fscanf(fp,'%f',1);
p=fscanf(fp,'%d',[1,n]);
fclose(fp);
