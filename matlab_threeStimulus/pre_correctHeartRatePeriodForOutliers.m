function [hpOut, thpOut] = pre_correctHeartRatePeriodForOutliers(rPeakTimes, t, y, handles)

    % Derived from HRVAS's " function [nibi,art]=correctEctopic(ibi,opt)"
    % http://sourceforge.net/projects/hrvas/

    [~, handles.flags] = init_DefaultSettings();
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'correctHRforOutliers.mat';
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

    % transpose inputs
    %t = t';
    %y = y';
    
    debugPlotOn = 0;

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
        else
            artPer=false(size(y,1),1);
        end

        if sum(strcmp('sd', handles.parameters.heart.ectopicOutlierMethod) == 1)
            artSD=locateOutliers(t,y,'sd',handles.parameters.heart.ectopicSdParam);            
        else
            artSD=false(size(y,1),1);
        end

        if sum(strcmp('median', handles.parameters.heart.ectopicOutlierMethod) == 1)
            artMed=locateOutliers(t,y,'median',handles.parameters.heart.ectopicMedianParam);            
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


    % plot(t, y, 'r', thpOut, hpOut, 'b')
