% findcode     -  codes to find. eg. {'tr1' 'tr3' 'tr6'}
% replacecode  -  codes to be used as a replacement. e.g. 'Target'
% deltatimems
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009
function EEG = replacecode(EEG, findcode, replacecode, deltatimems)
if nargin<1
      help replacecode
      return
end
fs = EEG.srate;
if nargin<4
else
      deltatimesam = round((deltatimems*fs)/1000); % delta in samples
end
%
% extracting codes from EEG.event.type
%
if ischar(EEG.event(1).type)
      currentcodes = {EEG.event.type}; %strings
else
      currentcodes = [EEG.event.type]; %numeric code
end
%
% search for findcodes
%
if ~iscell(findcode)
      if isnumeric(findcode)
            findcode = num2cell(findcode);
      else
            findcode = cellstr(findcode);
      end
end
if ischar(findcode{1}) && iscell(currentcodes)
      indxfindcode  = strmatch(findcode(1), currentcodes, 'exact');
elseif ~ischar(findcode{1}) && ~iscell(currentcodes) % numeric
      indxfindcode  = find(currentcodes==findcode{1});
elseif ischar(findcode{1}) && ~iscell(currentcodes)
      numt = str2num(findcode{1});
      if ~isempty(numt)
            indxfindcode  = find(currentcodes==numt);
      else
            msgboxText = 'You specified string(s) as event code, but your current events are numeric.';
            title = 'ERPLAB: input format error';
            errorfound(msgboxText, title);
            return
      end
elseif ~ischar(findcode{1}) && iscell(currentcodes)
      indxfindcode  = strmatch(num2str(findcode{1}), currentcodes, 'exact');
end
nf = length(findcode);
if nf==1
      [EEG.event(indxfindcode).type] = deal(replacecode);
else
      for i=1:length(indxfindcode)
            indx = indxfindcode(i);
            if indx>1 && (indx+nf-1)<=length(EEG.event)
                  if max(diff([EEG.event(indx:indx+nf-1).latency]))<=2 && diff([EEG.event(indx-1:indx).latency])>2
                        [EEG.event(indx).type] = deal(replacecode);
                        EEG.event(indx+1:indx+nf-1) = [];
                  end
            elseif indx==1
                  if max(diff([EEG.event(indx:indx+nf-1).latency]))<=2
                        [EEG.event(indx).type] = deal(replacecode);
                        EEG.event(indx+1:indx+nf-1) = [];
                  end
            end
      end
end
EEG = eeg_checkset( EEG );
