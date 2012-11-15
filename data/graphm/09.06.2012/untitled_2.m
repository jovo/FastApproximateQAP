for i = 1:137;
    
    A = a(i).G;
    B = a(i).H;

    I = eye(size(A,1));

    P = I(a(i).I,:);
    if ( ~isempty(P) )
        v = trace(B*P'*A*P);
        a(i).objFcns.I = v;
    else
        a(i).objFcns.I = nan;
    end
    
    P = I(a(i).U,:);
    if ( ~isempty(P) )
        v = trace(B*P'*A*P);
        a(i).objFcns.U = v;
    else
        a(i).objFcns.U = nan;
    end
    
    P = I(a(i).RANK,:);
    if ( ~isempty(P) )
        v = trace(B*P'*A*P);
        a(i).objFcns.RANK = v;
    else
        a(i).objFcns.RANK = nan;
    end
    
    P = I(a(i).QCV,:);
    if ( ~isempty(P) )
        v = trace(B*P'*A*P);
        a(i).objFcns.QCV = v;
    else
        a(i).objFcns.QCV = nan;
    end
    
    P = I(a(i).rand,:);
    if ( ~isempty(P) )
        v = trace(B*P'*A*P);
        a(i).objFcns.rand = v;
    else
        a(i).objFcns.rand = nan;
    end
    
    P = I(a(i).PATH,:);
    if ( ~isempty(P) )
        v = trace(B*P'*A*P);
        a(i).objFcns.PATH = v;
    else
        a(i).objFcns.PATH = nan;
    end
    
end