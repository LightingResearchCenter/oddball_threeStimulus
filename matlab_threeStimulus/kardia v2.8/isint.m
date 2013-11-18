%% ---- isint
function [d,ndx,varout] = isint(varargin)

    C = varargin;
    nIn = nargin;

    if nIn < 2
        d = C{:} == round(C{:});
    else
        d = zeros(1,nIn);
        for i = 1:nIn
            d(i) = isequal(C{i},round(C{i}));
        end
    end

    ndx = find(d);
    temp = [C{:}];
    varout = temp(ndx);
end

