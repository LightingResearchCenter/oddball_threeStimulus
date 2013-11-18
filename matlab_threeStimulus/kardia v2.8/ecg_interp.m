%% ---- ecg_interp
function yy=ecg_interp(t,y,tt,method)
    % ECG_INTERP interpolate HR series
    %   yy = ECG_INTERP(t,y,tt,method)
    %
    % Input arguments:
    %   t      - vector of HR instants
    %   y      - HR values
    %   tt     - vector of instants where HR values should be estimated
    %   method - 'constant', linear' or 'spline'
    %
    % Output arguments:
    %   yy     - HR values should at 'tt' instants
    %
    % Description:
    %   'constant' seems to be the most frequently used technique to calculate
    %   the averaged instantaneous heart period or rate in psychophysiology
    %   studies. It considers that HR is constant between two adjacent HR intervals.
    %   However, it seems to be the least effective estimate of the
    %   true espectrum. Preference should be done to 'spline' for spectral
    %   estimation.
    %
    % References:
    %   Guimaraes & Santos (1998)
    % -------------------------------------------------------------------------
    % Written by Mateus Joffily - NeuroII/UFRJ & CNC/CNRS

    if length(y) ~= length(t)
        disp ([mfilename ': "y" and "t" vectors must have the same length.']);
        yy=NaN;
        return
    end

    switch method
        case 'constant'
            yy=[];
            for i=1:numel(tt)
                k=find((tt(i)-t)>0);
                if isempty(k)
                    yy=[yy NaN];
                else
                    yy=[yy y(k(end))];
                end
            end

        case 'linear'
            yy = interp1(t,y,tt,'linear');

        case 'spline'
            yy=spline(t,y,tt);

        otherwise
            yy=NaN;
    end
end