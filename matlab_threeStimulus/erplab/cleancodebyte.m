function EEG = cleancodebyte(EEG,byte)
% When you get unexpected codes in EEGLAB, from Biosemi recordings, you will
% be able to use cleancodebyte.m in order to ignore either the less significative
% byte (1rst byte) or the most significative byte (2nd byte) from your datasets.
%
% xxxxxxxx xxxxxxxx
% 2nd byte 1rst byte
%
% Example: ignore (clean) the second byte (most significative byte)
% >> EEG = cleancodebyte(EEG,2);
%
% Author: Javier Lopez-Calderon & Johanna Kreither
% Center for Mind and Brain
% University of California, Davis
% April 2009

if ~ismember(byte,[1 2])
    disp('Error: "BYTE" HAVE TO BE EITHER 1 OR 2')
    return
end

type_binary = cellstr(dec2bin([EEG.event.type]));
nevents     = length(EEG.event);

for i=1:nevents
    currbin = type_binary{i};
    if length(currbin)==16
        if byte==1
            EEG.event(i).type = bin2dec(currbin(1:8)); % ignore first byte, use the 2nd byte
        elseif byte==2
            EEG.event(i).type = bin2dec(currbin(9:16)); % ignore second byte, use the 1rst byte
        end
    else
        fprintf('Event code %g (%g) is a one byte number', i, EEG.event(i).type);
    end
end

