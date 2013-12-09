% Author Javier Lopez-Calderon
function Seterplabfontsize

fs = erpworkingmemory('fontsizeGUI');
fu = erpworkingmemory('fontunitsGUI');

%
% Call GUI
%
answer = setvalueGUI({fs fu});

if isempty(answer)
      disp('User selected Cancel')
      return
end

fontsizeGUI  = answer{1};
fontunitsGUI = answer{2};
erpworkingmemory('fontsizeGUI', fontsizeGUI);
erpworkingmemory('fontunitsGUI', fontunitsGUI);