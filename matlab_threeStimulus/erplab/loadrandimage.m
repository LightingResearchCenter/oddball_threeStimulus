function [IconData IconCMap]= loadrandimage(fnames)
% Author: Javier
% 2008
%
% Just for fun!
try
        indxerrorfig = evalin('base', 'indxerrorfig');
        if isempty(indxerrorfig)
                indxerrorfig = 1:length(fnames);
        end
catch
        indxerrorfig = 1:length(fnames);
end
l = length(indxerrorfig);
ranindx = randperm(l);
[IconData IconCMap]= imread(fnames{indxerrorfig(ranindx(1))});
assignin('base', 'indxerrorfig', indxerrorfig(ranindx(2:end)));
return
