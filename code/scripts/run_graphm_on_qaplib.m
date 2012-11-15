% algNames = {'I','U','RANK','QCV','rand','PATH'};
algNames = {'PATH'};

for i=1:length(algNames)
    
    algName = algNames{i}
    
    [mn,s1,s3,s10]=runGraphmOnQAPs(algName,'all');
    
    nowstamp = strrep(strrep(datestr(now),' ','_'),':','.');
    
    save( [algName '_' nowstamp '.mat' ] );
    
end