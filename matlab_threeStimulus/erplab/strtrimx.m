% Remove leading and trailing white space from string, and remove any
% duplicated white space in between the strings (if any)  
% Javier Lopez-Calderon

function string = strtrimx(string)
string = strtrim(string);
string = regexprep(string, '\s+',' ');
