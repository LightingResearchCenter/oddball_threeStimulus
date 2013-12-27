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
    
        % ERP Waveforms (~400 MB each)
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_fullEpochs.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_epochs = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % ERP Waveforms (Stat, ~550 kB each)
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_statEpochs.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_statEpochs = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
    
        % Derived ERP Measures such as amplitude latency etc.
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_analyzed.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_analyzed = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % Time-Frequency derived (change the filenames when you run MAIN_ .. again)
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_waveletEpochs.mat')); % Specifies the type of files to be listed        
        %handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_waveletDerived.mat')); % Specifies the type of files to be listed     
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_waveletDerived = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % Time-Frequency, all epochs 
        %{
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_waveletEpochs.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_waveletEpochs = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        %}
        
        % EEG Fractal Analysis
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_fractalEEG.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_fractalEEG = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % "Extra measures", such as EOG and HEART data
        handles.batch.dirOutput     = dir(fullfile(handles.path.matFilesOut, '*0_extraAnalysis.mat')); % Specifies the type of files to be listed        
        handles.batch.fileNames     = {handles.batch.dirOutput.name}'; % Prints the filenames from the input folder into a variable                
        handles.batch.fileNameFields_extra = batch_separateFileNamesToFields(handles.batch.fileNames); % separate into fields
        
        % save finally the "raw filenames" if you wanna manually start
        % excluding data from a specific subject
        fileNames = handles.batch.fileNames;
        for i = 1 : length(fileNames)
            fileNames{i} = strrep(fileNames{i}, '_extraAnalysis.mat', '');
        end
        
    %% SPECIFY OUTLIER DATA (bad sessions, etc.)
    
        outlierFilenameList = {'ha_02_10.bdf'; 'sadasdha_02_10.bdf'};
        
    %% Compare conditions (session, intensity, subject) for 1D signals 
            
        % i.e. 1D signals such as the ERP components, Reaction time,
        % Band powers and IAF        
        batchCompare_componentMain(handles.batch.fileNameFields_analyzed, fileNames, outlierFilenameList, handles)
        
    %% Compare extra sensors
    
        batchCompare_extraSensors(handles.batch.fileNameFields_extra, handles)            
        
    %% Compare conditions for 2D signals (time-frequency)
    
        % batchCompare_timeFreq(handles.batch.fileNameFields_waveletDerived, handles)
        
    %% EPOCH COMPARISON
    
        % implement later
        % batchCompare_epochWaveforms(handles.batch.fileNameFields_epochs, outlierFilenameList, 'full', handles)
        
        batchCompare_epochWaveforms(handles.batch.fileNameFields_statEpochs, outlierFilenameList, 'stat', handles)
        
        
        
        
        
    