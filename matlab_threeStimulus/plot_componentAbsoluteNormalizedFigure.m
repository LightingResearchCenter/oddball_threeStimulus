% ABSOLUTE NORMALIZED
function plot_componentAbsoluteNormalizedFigure(fig, statsOut, normalizationTypes, erpBandType, erpComponent, erpFilterType, chsToPlot, erpTypes, fieldValue, subjects, handles)

    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    % normalizationTypes = {'rawWithNoNormalization'}
    rows = length(normalizationTypes);
    if strcmp(erpComponent, 'RT')
        cols = 1;
        set(fig, 'Position', [0.05*scrsz(3) 0.05*scrsz(4) 0.32*scrsz(3) 0.9*scrsz(4)])
    else
        cols = length(chsToPlot) * length(erpTypes);
        set(fig, 'Position', [0.05*scrsz(3) 0.05*scrsz(4) 0.92*scrsz(3) 0.9*scrsz(4)])
    end


    for normType = 1 : length(normalizationTypes)

        if strcmp(erpComponent, 'RT')
            index = normType;
            stim = 1; % target
            sp(index) = subplot(rows,cols,index);
            [p(index,:), styleHandles(index,:), yLims(index,:)] = plot_sessionMeanSubplot(statsOut.(normalizationTypes{normType}), ...
                        'Cz', normalizationTypes, erpBandType, erpTypes, rows, cols, normType, stim, [], index, fieldValue, erpComponent, erpFilterType, subjects, handles);
            drawnow
        else
            for stim = 1 : length(erpTypes)   
                for ch = 1 : length(chsToPlot)
                    index = ((normType-1) * cols) + ((stim-1) * length(chsToPlot)) + ch;
                    sp(index) = subplot(rows,cols,index);
                    % stimulusTypes{stim}
                    [p(index,:), styleHandles(index,:), yLims(index,:)] = plot_sessionMeanSubplot(statsOut.(normalizationTypes{normType}), ...
                        chsToPlot{ch}, normalizationTypes, erpBandType, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpFilterType, subjects, handles);
                    drawnow
                end
            end
        end

    end

    if strcmp(erpComponent, 'RT')
        set(sp(1), 'YLim', [200 800])
        set(sp(2:end), 'YLim', [-0.6 0.6])            

    elseif strcmp(fieldValue, 'peakMeanAmplit')

        if strcmp(erpComponent, 'P3') || strcmp(erpComponent, 'CNV')
            set(sp(1:cols), 'YLim', [0 30])
            set(sp(cols+1:(2*cols)), 'YLim', [-1 1])
            set(sp((2*cols)+1:end), 'YLim', [-1 1])
        elseif strcmp(erpComponent, 'N2') || strcmp(erpComponent, 'N1')
            set(sp(1:cols), 'YLim', [-10 10])
            set(sp(cols+1:(2*cols)), 'YLim', [-1 1])
            set(sp((2*cols)+1:end), 'YLim', [-1 1])
        end

    else

    end

    set(sp, 'XLim', [0.5 4.5], 'XTick', [1 2 3 4])

    try
        if handles.figureOut.ON == 1    
            drawnow
            dateStr = plot_getDateString(); % get current date as string
            %cd(path.outputFigures)            
            fileNameOut = ['plot_componentComparison_',  erpComponent, '_', erpFilterType, '_', fieldValue, '_', chsToPlot{1}, '-',  chsToPlot{2}, '_', erpBandType, '_', dateStr];
            export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            %cd(path.code)
        end
    catch err            
        if strcmp(err.identifier, '????')
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        else
            err
            err.identifier
        end
    end




