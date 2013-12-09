% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2011

function [WinRej2 ChanRej2 ] = discardshortsegments(WinRej, chanrej, shortsegsam, dwarning)

WinRej2= []; ChanRej2 = [];
if nargin<4
        dwarning = 1;
end
if dwarning
        fprintf('\nWARNING: Marked segments that are shorter than %g samples will unmarked.\n\n', shortsegsam);
end
nwin    = size(WinRej,1);
indxgood = [];
k=1;
for j=1:nwin
        widthw = WinRej(j,2) - WinRej(j,1);
        if widthw>shortsegsam
                indxgood(k) = j;
                k=k+1;
        end
end
WinRej2 = WinRej(indxgood,:);
if ~isempty(chanrej)
        ChanRej2= chanrej(indxgood,:);
end