function [f,myp,iter] = graphm(A,B,algName,IMAX,rootDir,qaplibp)
% On output:
%     f=sum(sum(A.*B(p,p))), where
%     p is the permutation found by FW after projecting the interior point
%         to the boundary.
%     iter is the number of iterations performed.

file1 = [rootDir '/' 'graphm/' 'temp/' 'A.dat'];
file2 = [rootDir '/' 'graphm/' 'temp/' 'B.dat'];

configFileIn = [rootDir '/' 'graphm/' 'temp/' 'graphm.config'];
expOutFile = [rootDir '/' 'graphm/' 'temp/' 'graphm.exp_out'];
verboseOutFile = [rootDir '/' 'graphm/' 'temp/' 'graphm.verbose'];

if strcmpi(algName,'PATH')
    A = max(max(A)) - A;
    B = max(max(B)) - B;
end

system(['rm ' expOutFile]);

writeAdjToGraphm(A,file1);
writeAdjToGraphm(B,file2);

cmdString = [rootDir '/' 'graphm/' 'bin/graphm' ' ' configFileIn ' "algo=' algName ' s;algo_init_sol=unif s;graph_1=' file1 ' s;graph_2=' file2 ' s;exp_out_file=' expOutFile  ' s;verbose_file=' verboseOutFile  ' s"'];

[status,result] = system( cmdString );

p = graphmReadPermutation(algName, expOutFile);

f = sum(sum(A.*B(p,p)));

myp = p;

iter = graphmReadIterations(algName, verboseOutFile);


end