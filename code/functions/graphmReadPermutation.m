function myp = graphmReadPermutation(algName, pathToExpOutFile)

text = fileread(pathToExpOutFile);

s = regexprep(text, ' \n', ',');

[~,~,~,matchstring,~,~,~] = regexp(s,[algName '(,[0-9]+)+']);

origStr = matchstring{1};
modifiedStr = strrep(origStr, [algName ','], '');
p = eval(['[' modifiedStr ']']);
    
myp = p;

end