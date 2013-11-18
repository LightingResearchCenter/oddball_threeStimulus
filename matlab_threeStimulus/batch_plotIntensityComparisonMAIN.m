function batch_plotIntensityComparisonMAIN(statsOut, statsPer, erpComponent, erpDataType, fieldValue, fileNameFields, stimulusType, chsToPlot, handles)

    %% DEBUG
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
                        noOfSession = length(statsOut.darkCondition.standard.dark.Cz.mean);
    
    
    % Plot
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color', 'w');
            
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
                set(sp(cols+1:(2*cols)), 'YLim', [-0.6 0.6])
                set(sp((2*cols)+1:end), 'YLim', [-0.6 0.6])
            elseif strcmp(erpComponent, 'N2') || strcmp(erpComponent, 'N1')
                set(sp(1:cols), 'YLim', [-10 10])
                set(sp(cols+1:(2*cols)), 'YLim', [-0.6 0.6])
                set(sp((2*cols)+1:end), 'YLim', [-0.6 0.6])
            end
        
        else
            
        end
        
        % yLims

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
    
    
    function [p, styleHandles, yLims] = plot_sessionMeanSubplot(statData, chToPlot, normalizationTypes, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpDataType, handles)

        i = 1;

        x = 1:1:4;
        
        % displace the x-axis a bit for better visualization
        x_displ = 0.1;
        x1 = x;
        x2 = x - x_displ;
        x3 = x + x_displ;
        
        % statData
        % statData.target
        % erpTypes{stim}
        % statData.(erpTypes{stim})
        if strcmp(erpComponent, 'RT')       
            if strcmp(normalizationTypes{normType}, 'absolute')
                secondsToMilliseconds = 1000;
            else
                secondsToMilliseconds = 1;
            end            
            try
                y1 = secondsToMilliseconds * statData.(erpTypes{stim}).dark.(chToPlot).mean;
            catch err
                err
                statData
                erpTypes{stim}
                statData.(erpTypes{stim})
                statData.(erpTypes{stim}).dark
                statData.(erpTypes{stim}).dark.(chToPlot)
                error('indexing incorrect')
            end
            y2 = secondsToMilliseconds * statData.(erpTypes{stim}).dim.(chToPlot).mean;
            y3 = secondsToMilliseconds * statData.(erpTypes{stim}).bright.(chToPlot).mean;
            
            % these errors are SDs of the means, thus for normalized values
            % you won't have by definition any error, add later the LE and
            % UE obtained in batch_normalizeComponents() using the
            % low-level subfunction batch_normalizeLowLevel()
            err1 = secondsToMilliseconds * statData.(erpTypes{stim}).dark.(chToPlot).SD;
            err2 = secondsToMilliseconds * statData.(erpTypes{stim}).dim.(chToPlot).SD;
            err3 = secondsToMilliseconds * statData.(erpTypes{stim}).bright.(chToPlot).SD;
        else
            y1 = statData.(erpTypes{stim}).dark.(chToPlot).mean;
            y2 = statData.(erpTypes{stim}).dim.(chToPlot).mean;
            y3 = statData.(erpTypes{stim}).bright.(chToPlot).mean;
            
            err1 = statData.(erpTypes{stim}).dark.(chToPlot).SD;
            err2 = statData.(erpTypes{stim}).dim.(chToPlot).SD;
            err3 = statData.(erpTypes{stim}).bright.(chToPlot).SD;
        end
        
        % PLOT
        hold on
        p(1) = errorbar(x1, y1, err1, 'ko');
        p(2) = errorbar(x2, y2, err2, 'r*');
        p(3) = errorbar(x3, y3, err3, 'ro');
        hold off
            
        % Annotations
        
            % title & ylabel
            if rem(index-1,cols) == 0 || index == 1 % for each row
                normString = normalizationTypes{normType};
                
                if strcmp('absolute', normString)
                    if strcmp(erpComponent, 'RT')
                        yLabelString = 'Latency [ms]';
                    else
                        yLabelString = '\muV';
                    end
                else
                    yLabelString = '\Delta';
                end                
            elseif (index - (floor(index/cols)*cols)) == 2
                normString = erpComponent;
                yLabelString = ' ';
            elseif (index - (floor(index/cols)*cols)) == 4
                normString = fieldValue;                
                yLabelString = ' ';
            elseif (index - (floor(index/cols)*cols)) == 0
                normString = erpDataType;                
                yLabelString = ' ';                
            else                
                normString = ' ';
                yLabelString = ' ';
            end

            if rem(index, length(chToPlot)) ~= 0
                stimString = erpTypes{stim};
            else
                stimString = ' ';
            end
            if strcmp(erpComponent, 'RT')
                titString = sprintf('%s', normString);    
            else
                titString = sprintf('%s\n%s\n%s', normString, stimString, chToPlot);    
            end
            styleHandles.tit = title(titString);
            
            % xlabel
            if (rows - 1)*cols < index
                xLabString = sprintf('%s', 'Session');
            else
                xLabString = ' ';
            end
            styleHandles.xLab = xlabel(xLabString);
        
            % ylabel
            styleHandles.yLab = ylabel(yLabelString);
            
            % LEGEND
            if strcmp(erpComponent, 'RT')
                if index == 1
                    leg = legend('dark', '10lux', '40lux', 3);
                        %legend('boxoff')
                        set(leg, 'Location', 'NorthWest')
                        set(leg, 'Position', [0.807142193308545 0.911998913043478 0.160780669144981 0.0654891304347826])
                else
                    leg = legend(' ', ' ', ' ', 3);
                        legend('hide')
                end
                
            else
                if (index - (floor(index/cols)*cols)) == 0
                    leg = legend('dark', '10lux', '40lux', 3);
                        %legend('boxoff')
                        set(leg, 'Location', 'NorthWest')
                        set(leg, 'Position', [0.9034 0.2538 0.050776 0.04918])
                else
                    leg = legend(' ', ' ', ' ', 3);
                        legend('hide')
                end
            end
            
            
        % Style
        set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)
        set(styleHandles.tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
            set(styleHandles.tit, 'HorizontalAlignment', 'center')
        set(styleHandles.xLab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
        set(styleHandles.yLab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
        set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
            set(leg, 'EdgeColor', [.4 .4 .4])
        
        % plots
        set(p(1),'MarkerFaceColor',[0 0 0],'Color',[0 0 0], 'Marker', 's');
        set(p(2),'MarkerFaceColor',[0 0.20 0.80],...
            'Marker','pentagram',...
            'Color',[0 0.20 0.80]);
        set(p(3),'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
        set(p, 'MarkerSize', handles.style.markerSize+1)

        
        
        
        % get y limits
        yLims = [min(min(([y1 y2 y3]))) max(max(([y1 y2 y3])))];
            
        