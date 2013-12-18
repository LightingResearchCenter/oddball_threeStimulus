function batchCompare_extraSensors(fileNameFields, handles)

    %% DEBUG
    debugMatFileName = 'tempComponentCompareExtra.mat';
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
        

    % Pull out the data   
    [heart, eye] = batch_pullOut_extraSensors(fileNameFields, handles);
    fractalEEG = [];
    % fractalEEG = .. use another subfunction as saved to other .mat files    
    
    % Do the stats for those component
    % statsPer = 'trials';
    statsPer = 'session'; % average one session        
    stimulusType = {'target'; 'distracter'; 'standard'};                
    [heart_Out, fractalEEG_Out, eye_Out] = batch_extraSensorStats(heart, fractalEEG, eye, statsPer, stimulusType, handles);
    
    % Plot the results
    plot_EXTRA_MAIN(heart_Out, fractalEEG_Out, eye_Out, handles)
    