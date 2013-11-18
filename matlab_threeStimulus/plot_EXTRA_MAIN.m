function plot_EXTRA_MAIN(heart_Out, fractalEEG_Out, eye_Out, handles)

    %% DEBUG
    debugMatFileName = 'tempExtraPlot.mat';
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
    
    
    fieldsOfInterest = {'HR_Mean'; 'NHF'; 'avIBI'; 'RMSSD'; 'DFA_alphaScaling'; 'MFDFA_mean_h'};
    %batch_plotAuxSubLoop(heart_Out, 'scalar', fieldsOfInterest, 'Heart', handles) 
    
    fieldsOfInterest = {'MFDFA_mean_h'; 'MFDFA_mean_Dh'; 'MFDFA_width_h';  'MFDFA_height_Dh'};
    batch_plotAuxSubLoop(fractalEEG_Out, 'scalar', fieldsOfInterest, 'fractalEEG', handles) 
    
    
    function batch_plotAuxSubLoop(extraStat, plotType, fieldsOfInterest, sensorType, handles)      
        
        scrsz = get(0,'ScreenSize'); % get screen size for plotting
        
        % init
        normFields = fieldnames(extraStat);
        conditions = fieldnames(extraStat.(normFields{1}));        
        dataTypes = fieldnames(extraStat.(normFields{1}).(conditions{1}).(plotType))
        
        cols = length(fieldsOfInterest);
        
        fig(1) = figure('Color', 'w', 'Name', sensorType);
        offset = 0.015;
            set(fig(1), 'Position', [(0.05+offset)*scrsz(3) (0.20-offset)*scrsz(4) 0.92*scrsz(3) 0.62*scrsz(4)])
            
        plottableFieldsFound = 0;
        
        for field = 1 : length(dataTypes)            
            
            if sum(strcmp(dataTypes{field}, fieldsOfInterest) == 1)            
                
                plottableFieldsFound = plottableFieldsFound + 1;
                
                for i = 1 : length(normFields)        

                    statParam = fieldnames(extraStat.(normFields{i}).(conditions{1}).(plotType).(dataTypes{field}));                    

                        for cond = 1 : length(conditions)                                            

                            % [param i cond field]

                            % in
                            structIn = extraStat.(normFields{i}).(conditions{cond}).(plotType).(dataTypes{field});

                            % out
                            extraStat.(normFields{i}).(conditions{cond}).(plotType).(dataTypes{field})
                            meanValue.(conditions{cond}) = extraStat.(normFields{i}).(conditions{cond}).(plotType).(dataTypes{field}).mean;
                            SD.(conditions{cond}) = extraStat.(normFields{i}).(conditions{cond}).(plotType).(dataTypes{field}).SD;
                            n.(conditions{cond}) = extraStat.(normFields{i}).(conditions{cond}).(plotType).(dataTypes{field}).n;

                            % debug                        
                            % debug = [meanValue.(conditions{cond}) SD.(conditions{cond}) n.(conditions{cond})]                  

                        end

                        rows = length(normFields);                    

                        % plot
                        index = (i-1)*cols + plottableFieldsFound;
                        sp(index) = subplot(rows, cols, index);
                        [p(index,:), styleHandles(index), yLims(index,:)] = plot_auxSubplot(meanValue, SD, n, index, rows, cols, handles, [], dataTypes{field}, normFields{i});
                        drawnow


                end

                try
                    if handles.figureOut.ON == 1    
                        drawnow
                        dateStr = plot_getDateString(); % get current date as string
                        %cd(path.outputFigures)            
                        fileNameOut = ['plot_extraSensors_',  sensorType, '_', dateStr];
                        disp([' ... saving figure to disk (', fileNameOut, '.png]'])
                        export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig(1))
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
            else
                % just don't plot that variable
                
            end
        end
        
    function [p, styleHandles, yLims] = plot_auxSubplot(meanValue, SD, n, index, rows, cols, handles, statParam, dataType, normField)
        
        x = 1:1:4;
      
        yLims = [1 2];
        
        % displace the x-axis a bit for better visualization
        x_displ = 0.1;
        x1 = x;
        x2 = x - x_displ;
        x3 = x + x_displ;
        
        % y
        y1 = meanValue.dark;
        y2 = meanValue.dim;
        y3 = meanValue.bright;

        % error
        err1 = SD.dark;
        err2 = SD.dim;
        err3 = SD.bright;
        
         % PLOT
        hold on
        p(1) = errorbar(x1, y1, err1, 'ko');
        p(2) = errorbar(x2, y2, err2, 'r*');
        p(3) = errorbar(x3, y3, err3, 'ro');
        hold off

        % Annotations
        
            % title & ylabel
            normString = statParam;
            if strcmp(normField, 'absolute')
                
                yLabelString = '?';
                
                if strcmp(statParam, 'IAF_amplitGravity')
                    yLabelString = 'Hz';
                end
                
                if strcmp(dataType, 'Amplit')
                    yLabelString = 'Power [mean, dB?]';
                elseif strcmp(dataType, 'PSD')
                    yLabelString = 'Power [integral, dB?]';
                elseif strcmp(dataType, 'ratio')
                    yLabelString = 'Ratio';
                end
                
            elseif strcmp(normField, 'darkCondition') || strcmp(normField, 'firstSession')
                yLabelString = '\Delta';
            else
                yLabelString = '???';
            end
                
            titString = sprintf('%s\n%s\n%s', normString, dataType, ' ');    
            styleHandles.tit = title(titString, 'Interpreter', 'none');
            
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
            if strcmp(statParam, 'RT')
                if index == 1
                    leg = legend('dark', '10lux', '40lux', 3);
                        %legend('boxoff')
                        set(leg, 'Location', 'NorthWest')
                else
                    leg = legend(' ', ' ', ' ', 3);
                        legend('hide')
                end
                
            else
                if (index - (floor(index/cols)*cols)) == 0
                    leg = legend('dark', '10lux', '40lux', 3);
                        %legend('boxoff')
                        set(leg, 'Location', 'NorthWest')
                        if strcmp(dataType, 'scalars')
                           set(leg, 'Position',[0.7570 0.27696 0.22944 0.08669])
                        else
                            set(leg, 'Position',[0.918 0.228 0.050776 0.06951])
                        end
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
        set(p, 'MarkerSize', handles.style.markerSize)
