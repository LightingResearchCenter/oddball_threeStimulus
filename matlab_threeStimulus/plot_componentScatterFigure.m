% STAT PLOT
function plot_componentStatFigure(fig, statsOut, statsTests, matricesSessionNorm, noOfSessions, normalizationTypes, erpBandType, erpComponent, erpFilterType, chsToPlot, chSelected, erpTypes, fieldValue, subjects, handles)

    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    if chSelected == 1
        set(fig, 'Position', [0.06*scrsz(3) 0.04*scrsz(4) 0.9*scrsz(3) 0.9*scrsz(4)])            
    elseif chSelected == 2
        set(fig, 'Position', [0.07*scrsz(3) 0.03*scrsz(4) 0.9*scrsz(3) 0.9*scrsz(4)])
    end        
    set(fig, 'Name', ['Scatter, ', chsToPlot{chSelected}])

    scatterTypes = {'absolute'; 'firstSession'};
    normalizationTypes = {'firstSession'};
        rows = length(normalizationTypes) + length(scatterTypes);
        cols = length(chSelected) * length(erpTypes);

    uniqueSubjects = unique(subjects);

    for normType = 1 : length(scatterTypes)
        for stim = 1 : length(erpTypes)   
            for ch = 1 : length(chSelected)
                index = ((normType-1) * cols) + ((stim-1) * length(chSelected) + ch);
                sp(index) = subplot(rows,cols,index);
                    [p(index,:,:), styleHandles(index,:), yLims(index,:)] = plot_sessionScatterSubplot(matricesSessionNorm.(scatterTypes{normType}), ...
                        chsToPlot{chSelected}, scatterTypes, erpBandType, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpFilterType, noOfSessions, uniqueSubjects, handles);
                    drawnow
            end
        end

    end               

    i1 = normType;
    for normType = 1 : length(normalizationTypes)
        j = i1 + normType;
        for stim = 1 : length(erpTypes)   
            for ch = 1 : length(chSelected)
                index = ((j-1) * cols) + ((stim-1) * length(chSelected)) + ch;
                sp(index) = subplot(rows,cols,index);
                % stimulusTypes{stim}
                    [p2(normType,:), styleHandles(index,:), yLims(index,:)] = plot_sessionMeanSubplot(statsOut.(normalizationTypes{normType}), ...
                        chsToPlot{chSelected}, normalizationTypes, erpBandType, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpFilterType, uniqueSubjects, handles);
                    drawnow
            end
        end

    end

    if strcmp(fieldValue, 'peakMeanAmplit')

        if strcmp(erpComponent, 'P3') || strcmp(erpComponent, 'CNV')
            %set(sp(1:cols), 'YLim', [-10 30])
            set(sp(1:cols), 'YLim', [min(min(yLims(1:cols,:))) max(max(yLims(1:cols,:)))])
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
            fileNameOut = ['plot_componentScatter_',  erpComponent, '_', erpFilterType, '_', fieldValue, '_', chsToPlot{chSelected}, '_', erpBandType, '_', dateStr];
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

