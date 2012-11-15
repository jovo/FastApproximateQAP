for i = 1:822
    
    textdata(i,1) ; % algorithm
    textdata(i,2) ; % dataset
    data(i) ; % time
    
    for j = 1:137
        if ( strcmp(textdata(i,2), a(j).name ) )
            if ( strcmp(textdata(i,1), 'I') )
                a(j).times.I = data(i);
            elseif ( strcmp(textdata(i,1), 'U') )
                a(j).times.U = data(i);
            elseif ( strcmp(textdata(i,1), 'RANK') )
                a(j).times.RANK = data(i);
            elseif ( strcmp(textdata(i,1), 'QCV') )
                a(j).times.QCV = data(i);
            elseif ( strcmp(textdata(i,1), 'rand') )
                a(j).times.rand = data(i);
            elseif ( strcmp(textdata(i,1), 'PATH') )
                a(j).times.PATH = data(i);
            end
        end
    end 
end