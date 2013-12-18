function [epochsOut, debugInfoOut] = batch_pullOut_timeDomainEpochs(fileNameFields, outlierFilenameList, erpType, erpFilterType, epochType, handles)

    warning on
    
    %% DEBUG
    debugMatFileName = 'tempPullOutEPOCHS.mat';
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
    
    % check input
    numberOfFilesFound = length(fileNameFields);
    filesPerSubject = 12;
    numberOfSubjectsDone = numberOfFilesFound/filesPerSubject;
    disp([num2str(numberOfFilesFound), ' files found, thus most likely ', num2str(numberOfSubjectsDone), ' subjects done?'])
    
    if rem(numberOfSubjectsDone, 1) % if non-integer
       error('Lazy programmer have not managed yet the situation where you do not all the sessions done for each subject')       
    end
    
    % outlier files
    outlierFilesFoundFromInput = false(length(outlierFilenameList), 1);
    
    %% Go through the files
    for i = 1 : numberOfFilesFound                       
        
        % check intensity
        if fileNameFields{i}.intensity == 0
            intensity = 'dark';
        elseif fileNameFields{i}.intensity == 10
            intensity = 'dim';
        elseif fileNameFields{i}.intensity == 40
            intensity = 'bright';
        else
            intensity
            a = fileNameFields{i}.intensity
            error('Error in the intensity field!')
        end

        % check session
        if fileNameFields{i}.session == 1
            session = 'session1';
        elseif fileNameFields{i}.session == 2
            session = 'session2';
        elseif fileNameFields{i}.session == 3
            session = 'session3';
        elseif fileNameFields{i}.session == 4
            session = 'session4';
        else
            session
            error('Error in the session field!')
        end          
        
        % load the actual data
        disp(['    ... load the data from file: "', fileNameFields{i}.fileName, '" (', num2str(i), ' / ', num2str(numberOfFilesFound), ')'])
        dataIn = load(fullfile(handles.path.matFilesOut, fileNameFields{i}.fileName));
        
        % check if the input has been manually marked to be an outlier
        extraString = '_fullEpochs.mat';
        outlierListTrue = strcmp(strrep(fileNameFields{i}.fileName, extraString, '.bdf'), outlierFilenameList);
        thisFileMarkedAsOutlier = logical(sum(outlierListTrue));
        if thisFileMarkedAsOutlier == 1
            disp(['        - ', fileNameFields{i}.fileName, ' was marked to be an outlier'])
        end
        
            % check also that all the files marked as outliers are
            % found from the input files, i.e. if you have a typo in
            % your outlier files, it will be never found and corrected,
            % while you might think that it is been corrected
            outlierFilesFoundFromInput(outlierListTrue) = 1;
        
        % subject
        subject = fileNameFields{i}.subject;
        subjects{i} = subject;
        
        % should return target, distr, std
        responseTypes = fieldnames(dataIn.epochs.(erpFilterType).(erpType));
        
        % go through the different response types, and average the trials
        % together already here to speed up the followign processing
        for resp = 1 : length(responseTypes)
            
            % process using a subfunction
            [epochsProcessed, debug] = batch_timeDomEpochsAverage(dataIn.epochs.(erpFilterType).(erpType).(responseTypes{resp}), thisFileMarkedAsOutlier, handles.parameters, handles);
            
            % assign to output
            epochsOut.(intensity).(session).(responseTypes{resp}).(erpType).(subject).(erpFilterType) = epochsProcessed;
                    
            % assign info on trigger variability
            debugInfoOut.(subject).(intensity).(session) = debug;
            
        end
        
        % hh
        % dads
        parameters = handles.parameters;
        fileMatOut = strrep(fileNameFields{i}.fileName, '_fullEpochs.mat', '_statEpochs.mat');        
        disp(['      .. Stats of ERP Epochs [', fullfile(handles.path.matFilesOut, fileMatOut), ']'])
        save(fullfile(handles.path.matFilesOut, fileMatOut), 'epochsProcessed', 'debug', 'parameters', 'handles')
        
        
    end
    
    
    %% Check if all the outlier files have been found            
    notFoundLinearIndices = find(outlierFilesFoundFromInput == 0);
    for i = 1 : length(notFoundLinearIndices)   
        warning('Not all the outlier files were found!')
        fprintf('     ')
        for j = 1 : i
           fprintf(' '); % add white space
        end        
        fprintf(['"', outlierFilenameList{notFoundLinearIndices(i)}, '" not found from the input files, is this correct?\n'])        
    end
    fprintf('\n')
        
        
    