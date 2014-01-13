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
    
        % get your computer name
        [ret, compName] = system('hostname'); 
            %compName = lower(compName);
            compName = strtrim(compName); % remove whitespace
            % can be used to define custom paths for each user developing/using this
            % code simultaneously

            % you need to define the home folder for the INPUT DATA,
            % 'homeFolder', and then the output folder matFilesMain
            
        % main paths 
        % change these only if you move the data somewhere else                
        if strcmp(compName, 'Ubuntu64')
            handles.path.homeFolder = '/home/petteri/CopyShare/'; % Unix home folder, for EEG Data            
        elseif strcmp(compName, 'someoneElse')
            
        elseif strcmp(compName, 'someoneElse2')
            
        else
            error('Unknown computer name, define your paths!')
        end
        
        
        handles.path.codeFolder = mfilename('fullpath'); % Setting the path for the code
        handles.path.codeFolder = strrep(handles.path.codeFolder, 'init_DefaultSettings',''); % Removing the filename from the path

            handles.path.dataFolder = fullfile(handles.path.homeFolder, 'EEG-threeStimOddball');            
                handles.path.dataKSS = fullfile(handles.path.dataFolder, 'KSS');
            
            handles.path.figuresOut = fullfile(handles.path.codeFolder, 'figuresOut');
            handles.path.textOut = fullfile(handles.path.codeFolder, 'textFilesOut');
            
            handles.path.debugOut = fullfile(handles.path.figuresOut, 'debug');
            handles.path.debugHeartOut = fullfile(handles.path.debugOut, 'heart');
            handles.path.debugPreprocessing = fullfile(handles.path.debugOut, 'preProcessing');
            handles.path.debugFractal = fullfile(handles.path.debugOut, 'fractalEEG');
            handles.path.debugFASTER = fullfile(handles.path.debugOut, 'FASTER');
            
            % OUT
            if strcmp(compName, 'Ubuntu64')
                handles.path.matFilesMain = '/home/petteri/EEG-threeStim/';
            elseif strcmp(compName, 'someoneElse')
                handles.path.matFilesMain = '????????????????';
            else
                error('Unknown computer name, define your paths!')
            end
            
            % check that the folder actually exists
            try
                currDir = pwd;
                cd(handles.path.matFilesMain)
                cd(currDir)
            catch err
                err
                error('OUTPUT Directory does not exist. Note that the output .MAT files are stored outside the Cloud/Github path!')
                % in other words change the handles.path.matFilesMain to a
                % folder that exists
            end
            handles.path.matFilesOut = fullfile(handles.path.matFilesMain, 'matOut');
            handles.path.matFilesInput = fullfile(handles.path.matFilesMain, 'inputMat');
            handles.path.debugMATs = fullfile(handles.path.matFilesMain, 'debugMATs');
            
            % add try/catch for folders, first time use!

    %% PLOT STYLING
    
        handles.style.scrsz = get(0,'ScreenSize'); % get screen size for plotting
        set(0,'DefaultFigureColor','w')
    
        handles.style.fontName = 'Latin Modern Roman';
        handles.style.fontSizeBase = 10;     
        handles.style.markerSize = 6;  
        handles.style.markerFaceColor = [0 0.4 1];
        handles.style.fontGrey = [0.2 0.2 0.2];
        handles.style.lineGrey = [0.4 0.4 0.4];
        handles.style.markerEdgeColor = 'none';
        handles.style.ERP_yLimits = [-10 15];
        handles.style.RT_limits = [-200 800];
           
        % settings when auto-saving figures, see exportfig.m for more details        
        handles.figureOut.ON                = 1;
        handles.figureOut.debugON           = 1; % saving the under-the-hood plots
        handles.figureOut.resolution        = '-r150';  
        handles.figureOut.format            = 'png';        
        handles.figureOut.antialiasLevel    = '-a1';