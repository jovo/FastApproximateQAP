for ii = 1:1000
    p = randperm(279);
    dlmwrite(sprintf('perm_%d.txt',ii),p,' ');
    dlmwrite(sprintf('Agap_%d.txt',ii),full(Agap(p,p)),' ');
    dlmwrite(sprintf('Achem_%d.txt',ii),full(Achem(p,p)),' ');
end
