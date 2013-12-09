% Author: Javier Lopez-Calderon
function [MPD com] = getardetection(ERPLAB, cw)
com = '';
if nargin<2
      cw = 0; % do not comment on command window
end
try
      if iseegstruct(ERPLAB)
            MPD = pop_summary_AR_eeg_detection(ERPLAB, 'none');
      elseif iserpstruct(ERPLAB)
            MPD = []; % pending
      else
            MPD = [];
      end
catch
      MPD = [];
end
if cw==1
      if isempty(MPD)
            fprintf('\nCannot get the information you requested. Sorry\n\n');
      else
            fprintf('\nMean artifact detection ratio : %.1f\n\n', MPD);
      end
end
com = sprintf('MPD = getardetection(%s);', inputname(1));
return
