function batch_plotIntensityComparisonMAIN(statsOut, matricesSessionNorm, statsPer, erpComponent, erpDataType, fieldValue, fileNameFields, stimulusType, chsToPlot, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempBatchPlot.mat';
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
    end
    
    %{
    whos
    statsOut
    statsOut.darkCondition
    statsOut.darkCondition.target
    statsOut.darkCondition.target.dark
    statsOut.darkCondition.target.dark.Fz
    statsOut.darkCondition.target.dark.Fz.mean
    %}
    
    % Data "hidden" in structure        
    normalizationTypes = fieldnames(statsOut);
        erpTypes = fieldnames(statsOut.darkCondition);
            conditionTypes = fieldnames(statsOut.darkCondition.standard);
                chNames = fieldnames(statsOut.darkCondition.standard.dark);
                    statFields = fieldnames(statsOut.darkCondition.standard.dark.Cz);
                        noOfSessions = length(statsOut.darkCondition.standard.dark.Cz.mean);
    
    
    %% PLOT 1: Absolute / Normalization

        
        fig1 = figure('Color', 'w');
        plot_componentAbsoluteNormalizedFigure(fig1, statsOut, normalizationTypes, erpComponent, erpDataType, chsToPlot, erpTypes, fieldValue, handles)

    
    %% PLOT 2:        

        fig2 = figure('Color', 'w');
        chSelected = 1; % too clogged if you plot 2 channels to same figure
        plot_componentScatterFigure(fig2, statsOut, matricesSessionNorm, noOfSessions, normalizationTypes, erpComponent, erpDataType, chsToPlot, chSelected,  erpTypes, fieldValue, handles)

        fig3 = figure('Color', 'w');
        chSelected = 2; % too clogged if you plot 2 channels to same figure
        plot_componentScatterFigure(fig3, statsOut, matricesSessionNorm, noOfSessions, normalizationTypes, erpComponent, erpDataType, chsToPlot, chSelected, erpTypes, fieldValue, handles)

    
%% SUBFUNCTIONS

    function plot_componentScatterFigure(fig, statsOut, matricesSessionNorm, noOfSessions, normalizationTypes, erpComponent, erpDataType, chsToPlot, chSelected, erpTypes, fieldValue, handles)

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


        for normType = 1 : length(scatterTypes)
            for stim = 1 : length(erpTypes)   
                for ch = 1 : length(chSelected)
                    index = ((normType-1) * cols) + ((stim-1) * length(chSelected) + ch);
                    sp(index) = subplot(rows,cols,index);
                    [p(index,:,:), styleHandles(index,:), yLims(index,:)] = plot_sessionScatterSubplot(matricesSessionNorm.(scatterTypes{normType}), ...
                        chsToPlot{chSelected}, scatterTypes, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpDataType, noOfSessions, handles);
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
                        chsToPlot{chSelected}, normalizationTypes, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpDataType, handles);
                    drawnow
                end
            end

        end

        if strcmp(fieldValue, 'peakMeanAmplit')

            if strcmp(erpComponent, 'P3') || strcmp(erpComponent, 'CNV')
                set(sp(1:cols), 'YLim', [-10 30])
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
                fileNameOut = ['plot_componentScatter_',  erpComponent, '_', erpDataType, '_', fieldValue, '_', chsToPlot{chSelected}, '_', dateStr];
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

    
    
    function plot_componentAbsoluteNormalizedFigure(fig, statsOut, normalizationTypes, erpComponent, erpDataType, chsToPlot, erpTypes, fieldValue, handles)

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
                            'Cz', normalizationTypes, erpTypes, rows, cols, normType, stim, [], index, fieldValue, erpComponent, erpDataType, handles);
                drawnow
            else
                for stim = 1 : length(erpTypes)   
                    for ch = 1 : length(chsToPlot)
                        index = ((normType-1) * cols) + ((stim-1) * length(chsToPlot)) + ch;
                        sp(index) = subplot(rows,cols,index);
                        % stimulusTypes{stim}
                        [p(index,:), styleHandles(index,:), yLims(index,:)] = plot_sessionMeanSubplot(statsOut.(normalizationTypes{normType}), ...
                            chsToPlot{ch}, normalizationTypes, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpDataType, handles);
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
                fileNameOut = ['plot_componentComparison_',  erpComponent, '_', erpDataType, '_', fieldValue, '_', chsToPlot{1}, '-',  chsToPlot{2}, '_', dateStr];
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

    
    

