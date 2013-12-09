% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2011

function [WinRej2 ChanRej2 ] = joinclosesegments(WinRej, chanrej, shortisisam)

WinRej2= []; ChanRej2 = [];
fprintf('\nWARNING: Marked segments that are closer than %g samples will be join together.\n\n', shortisisam);
chanrej = uint8(chanrej);
a = WinRej(1,1);
b = WinRej(1,2);
m = 1;
working = 0;
chrej2 = uint8(zeros(1,size(chanrej,2)));
nwin = size(WinRej,1);
for j=2:nwin
    isi = WinRej(j,1) - WinRej(j-1,2);
    if isi<shortisisam
        b = WinRej(j,2);
        chrej2 = bitor(chrej2, bitor(chanrej(j,:),chanrej(j-1,:)));
        working = 1;
        if j==nwin
            WinRej2(m,:)  = [a b];
            ChanRej2(m,:) = chrej2;
        end
    else
        if working==1
            WinRej2(m,:)  = [a b];
            ChanRej2(m,:) = chrej2;
            %                     a = WinRej(j,1);
            working = 0;
        else
            WinRej2(m,:)  = [a b];
            ChanRej2(m,:) = chanrej(j-1);
            %                     a = WinRej(j,1);
            %                     b = WinRej(j,2);
        end
        a = WinRej(j,1);
        b = WinRej(j,2);
        chrej2 = uint8(zeros(1,size(chanrej,2)));
        m = m + 1;
    end
end
ChanRej2  = double(ChanRej2);
