for i=1:137
    FAQ.times(i)=A(i).FAQ.t;
    FAQ.obj(i)=A(i).FAQ.obj;
    
    PATH.times(i)=A(i).PATH.t;
    PATH.obj(i)=A(i).PATH.obj;
    Nvertices(i)=length(A(i).FAQ.p);
end

