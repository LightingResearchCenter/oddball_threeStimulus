function [handles, flags] = init_DefaultSettings()
           
    handles.style.scrsz = get(0,'ScreenSize'); % get screen size for plotting
    handles.style.scrsz = [1 1 1680 1050];
        
    % General debug switch
    handles.flags.showDebugMessages = 0;
    handles.flags.showDebugPlots = 0;
    handles.flags.saveDebugMATs = 1; % takes quite a lot space from HDD, but easier to develop, 
                                     % no need to be ON when just using this
                                     % script and not adding anything new
                                     % Takes maybe roughly 25% more time to
                                     % process everything with this option
                                     % ON (depends on the speed of your HDD
                                     % though)

    flags = handles.flags;

    %% PATHS

        % main paths 
        % change these only if you move the data somewhere else
        handles.path.homeFolder = '/home/petteri/CopyShare/'; % Unix home folder, for EEG Data
        handles.path.codeFolder = mfilename('fullpath'); % Setting the path for the code
        handles.path.codeFolder = strrep(handles.path.codeFolder, 'init_DefaultSettings',''); % Removing the filename from the path

            % derived pathnames
            % NO NEED to touch these unless you know what you are really doing
            handles.path.dataFolder = fullfile(handles.path.homeFolder, 'EEG-threeStimOddball');
            handles.path.debugMATs = fullfile(handles.path.dataFolder, 'debugMATs');
            
            handles.path.figuresOut = fullfile(handles.path.codeFolder, 'figuresOut');
            
            handles.path.debugOut = fullfile(handles.path.codeFolder, 'debugOut');
            handles.path.debugHeartOut = fullfile(handles.path.debugOut, 'HRV');
            handles.path.debugPreprocessing = fullfile(handles.path.debugOut, 'preProcessing');
            
            handles.path.matFilesOut = '/home/petteri/EEG-threeStim/matOut/';

    %% PLOT STYLING
    
        handles.style.scrsz = get(0,'ScreenSize'); % get screen size for plotting
        set(0,'DefaultFigureColor','w')
    
        handles.style.fontName = 'Latin Modern Roman';
        handles.style.fontSizeBase = 10;     
        handles.style.markerSize = 6;  
        handles.style.markerFaceColor = [0 0.4 1];
        handles.style.fontGrey = [0.2 0.2 0.2];
        handles.style.markerEdgeColor = 'none';
        handles.style.ERP_yLimits = [-10 15];
        handles.style.RT_limits = [-200 800];
           
        % settings when auto-saving figures, see exportfig.m for more details        
        handles.figureOut.ON                = 1;
        handles.figureOut.debugON           = 1; % saving the under-the-hood plots
        handles.figureOut.resolution        = '-r150';  
        handles.figureOut.format            = 'png';        
        handles.figureOut.antialiasLevel    = '-a1';