function debug_plotTF(scales, timep, freqs, realCoefs, imagCoefs, T, F, handles)

    %% DEBUG
    debugMatFileName = 'debugTFPlot.mat';
    if nargin == 0
        load('debugPath.mat')
        load(fullfile(path.debugMATs, debugMatFileName))
        close all
    else
        if handles.flags.saveDebugMATs == 1
            path = handles.path;
            save('debugPath.mat', 'path')
            save(fullfile(path.debugMATs, debugMatFileName))            
        end
    end
    
    whos
    
    subplot(2,1,1)
    contourf(T, F, squeeze(realCoefs(2,:,:)), 64, 'EdgeColor', 'none');
    xlabel('Time [ms]'); ylabel('Frequency [Hz]');
    colorbar
    
    subplot(2,1,2)
    contourf(T, F, squeeze(imagCoefs(2,:,:)), 64, 'EdgeColor', 'none');
    xlabel('Time [ms]'); ylabel('Frequency [Hz]');
    colorbar