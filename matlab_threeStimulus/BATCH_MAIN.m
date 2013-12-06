% MAIN FUNCTION for Learning Oddball EEG ANALYSIS
function BATCH_MAIN()

    tic

    % Petteri Teikari, petteri.teikari@gmail.com, 2013
    % Lighting Research Center, Rensselaer Polytechnic Institute, Troy, NY
    close all
    clear all
    
    %% General Settings
    
        % i.e. like where the folders are, fonts to be used, colors, etc.
        [handles, handles.flags] = init_DefaultSettings(); % use a subfunction        

    %% Parameters for ANALYSIS
    
        % i.e. like artifact rejection thresholds, filter cutoffs,
        % numerical parameters for EEG analysis
        handles.parameters = init_DefaultParameters(handles); % use a subfunction        
    
    %% Get directory listing of the processed files
    
        % For 1D signals
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_analyzed.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % For 1D signals
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_extraAnalysis.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_extra = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % For Epochs
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_fullEpochs.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_epochs = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
    %% Compare conditions (session, intensity, subject) for 1D signals 
            
        % i.e. 1D signals such as the ERP components, Reaction time,
        % Band powers and IAF        
        % batchCompare_componentMain(handles.batch.fileNameFields, handles)
        
    %% Compare extra sensors
    
        batchCompare_extraSensors(handles.batch.fileNameFields_extra, handles)
            
        
    %% Compare conditions for 2D signals (time-frequency)
    
        batchCompare_timeFreq(handles.batch.fileNameFields_extra, handles)
        
    %% EPOCH COMPARISON
    
        % implement later
        batchCompare_epochWaveforms(handles.batch.fileNameFields_epochs, handles)
        
        
        
        
        
    