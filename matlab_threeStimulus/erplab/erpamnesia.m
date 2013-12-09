%
% Author: Javier Lopez-Calderon
%
function erpamnesia

p = which('eegplugin_erplab');
p = p(1:findstr(p,'eegplugin_erplab.m')-1);
recycle on;
tempname = fullfile(p, 'memoryerp.erpm');
delete(tempname)
recycle off;
fprintf('\nERPLAB''s memory was wiped out!\n\n');