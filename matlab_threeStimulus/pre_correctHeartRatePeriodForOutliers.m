function [hpOut, thpOut] = pre_correctHeartRatePeriodForOutliers(rPeakTimes, t, y, callFrom, handles)

    % Derived from HRVAS's " function [nibi,art]=correctEctopic(ibi,opt)"
    % http://sourceforge.net/projects/hrvas/

    [~, handles.flags] = init_DefaultSettings();
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'correctHRforOutliers.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
            debugPlotOn = 1;
        else
            if handles.flags.saveDebugMATs == 1
                path = handles.path;
                save('debugPath.mat', 'path')
                save(fullfile(path.debugMATs, debugMatFileName))            
            end
            debugPlotOn = 1;
        end  
    end

    % transpose inputs, if needed, the vector need to be a column vector
    % (i.e. one column and x number of rows)
    if size(t,2) > size(t,1)
        t = t';        
    end
    
    if size(y,2) > size(y,1)
        y = y';        
    end    
    
    %% ACTUAL ALGORITHM
        
        if strcmp(handles.parameters.heart.outlierReplaceMethod, 'mean')
            handles.parameters.heart.outlierReplaceInput = 9; % 9 is default value for median from HRVAS
        elseif strcmp(handles.parameters.heart.outlierReplaceMethod, 'median')
            handles.parameters.heart.outlierReplaceInput = 5; % 5 is default value for median from HRVAS
        else
            handles.parameters.heart.outlierReplaceInput = NaN;
        end


        %% Locate ectopic

            if sum(strcmp('percent', handles.parameters.heart.ectopicOutlierMethod) == 1)
                artPer=locateOutliers(t,y,'percent',handles.parameters.heart.ectopicPercentParam);    
                paramDebug = handles.parameters.heart.ectopicPercentParam;
            else
                artPer=false(size(y,1),1);
            end

            if sum(strcmp('sd', handles.parameters.heart.ectopicOutlierMethod) == 1)
                artSD=locateOutliers(t,y,'sd',handles.parameters.heart.ectopicSdParam);            
                paramDebug = handles.parameters.heart.ectopicSdParam;
            else
                artSD=false(size(y,1),1);
            end

            if sum(strcmp('median', handles.parameters.heart.ectopicOutlierMethod) == 1)
                artMed=locateOutliers(t,y,'median',handles.parameters.heart.ectopicMedianParam);            
                paramDebug = handles.parameters.heart.ectopicMedianParam;
            else
                artMed=false(size(y,1),1);
            end

            art = artPer | artSD | artMed; %combine all logical arrays


        %% Replace ectopic

             switch handles.parameters.heart.outlierReplaceMethod
                case 'mean'
                    [hpOut, thpOut]=replaceOutliers(t,y,art,'mean',handles.parameters.heart.outlierReplaceInput);
                case 'median'
                    [hpOut, thpOut]=replaceOutliers(t,y,art,'median',handles.parameters.heart.outlierReplaceInput);
                case 'spline'
                    [hpOut, thpOut]=replaceOutliers(t,y,art,'cubic');
                case 'remove'
                    try
                        [hpOut, thpOut]=replaceOutliers(t,y,art,'remove');            
                    catch err
                        disp('           .. Had to transpose the t and y for artifact removal')
                        [hpOut, thpOut]=replaceOutliers(t',y',art,'remove');            
                    end

                otherwise %none
                    error(['? (typo in your replace method)']) 
             end
             
             
     %% DEBUG    
        
        if debugPlotOn == 1

            fig = figure('Color', 'w', 'Name', 'R / RR Outlier Debug');

                % do all the methods for comparison
                artPer = locateOutliers(t,y,'percent',handles.parameters.heart.ectopicPercentParam);
                artSD = locateOutliers(t,y,'sd',handles.parameters.heart.ectopicSdParam);
                artMed = locateOutliers(t,y,'median',handles.parameters.heart.ectopicMedianParam);

                m = mean(y); m = m /2;

                % plot input
                plot(t, y, 'k')

                % plot artifacts
                hold on
                plot(t, m*artPer, 'r', t, m*artSD, 'g', t, m*artMed, 'b')
                hold off
                lab(1) = xlabel('Time [s]');
                if strcmp(callFrom, 'rPeak')
                    lab(2) = ylabel('Amplitude / Outlier Logical');
                    titStr = sprintf('%s\n%s', 'Outlier removal for R peak amplitudes', ['method used = "', handles.parameters.heart.ectopicOutlierMethod, '", with param = ', num2str(paramDebug)]);
                elseif strcmp(callFrom, 'rInterval')
                    lab(2) = ylabel('Interval / Outlier Logical');
                    titStr = sprintf('%s\n%s', 'Outlier removal for RR timeseries', ['method used = "', handles.parameters.heart.ectopicOutlierMethod, '", with param = ', num2str(paramDebug)]);
                else
                    warning([callFrom, ' is in incorrect input string, needed only for plotting though'])
                end
                
                tit = title(titStr);

                leg = legend(callFrom, ['Outlier Percent, n = ', num2str(sum(artPer == 1))], ...
                                            ['Outlier SD, n = ', num2str(sum(artSD == 1))], ...
                                            ['Outlier Median, n = ', num2str(sum(artMed == 1))]);

                set(leg, 'Location', 'NorthOutside')
                legend('boxoff')
                
                set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
                set(gca, 'XLim', [min(t) max(t)]) 
                set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
                set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')  
                set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
                
                % AUTO-SAVE FIGURE
                try
                    if handles.figureOut.debugON == 1 
                        drawnow
                        dateStr = plot_getDateString(); % get current date as string
                        %cd(path.outputFigures)            
                        fileNameOut = [callFrom, 'OutlierRemoval_',strrep(handles.inputFile, '.bdf', ''), '_', dateStr];
                        disp(['         ... saving figure to disk (', fileNameOut, '.png]'])
                        export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig(1))
                        %cd(path.code)
                    end
                catch err
                    err
                end     

        end
