function batch_writeDataToDiskForStats(inputMatrix, statsOut, stat_results, stat_assumptions, normFieldName, erpBandType, erpComponent, erpFilterType, ...
                                        statField, stimType, chName, sessions, conditions, subjects, dataType, parameters, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatWriteToDisk.mat';
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
    
    [noOfConditions, noOfSessions, noOfSubjects] = size(inputMatrix);
    uniqueSubjects = unique(subjects);
    
    % define filename    
    dateStr = plot_getDateString(); % get current date as string    
    if strcmp(dataType, 'component')
        fileNameOut = [dataType, '_', erpFilterType, '_', erpComponent, '_', 'bandpassFilter-', erpBandType, '_', stimType, '_', normFieldName, '_', chName, '_', statField, '_', dateStr, '.txt'];
    elseif strcmp(dataType, 'AUX')
        fileNameOut = [dataType, '_', erpFilterType, '_', erpComponent, '_', 'bandpassFilter-', erpBandType, '_', stimType, '_', normFieldName, '_', chName, '_', statField, '_', dateStr, '.txt'];
    else
        error(['Your datatype is: ', dataType, '. There is nothing coded for this dataType, or is this a typo?'])
    end
    fileOutWithPath = fullfile(handles.path.textOut, fileNameOut);
    
    %% HEADER   
            
        % define header rows
        if strcmp(dataType, 'component')
            i = 1;
            headerRow{i} = ['Numerical data of the "', dataType, '"']; i = i + 1;
            headerRow{i} = ['===================================================']; i = i + 1;
            headerRow{i} = ['ERP filtered using "', erpFilterType, '" filtering']; i = i + 1;
            bandPassFilterFieldLo = ['bandPass_', erpBandType, '_loFreq']; i = i + 1;
            bandPassFilterFieldHi = ['bandPass_', erpBandType, '_hiFreq']; i = i + 1;
            headerRow{i} = ['The bandpass type is "', erpBandType, '" (', num2str(parameters.filter.(bandPassFilterFieldLo)), '-', num2str(parameters.filter.(bandPassFilterFieldHi)), ' Hz)']; i = i + 1;
            headerRow{i} = ['The audio tone type is: "', stimType, '"']; i = i + 1;
            headerRow{i} = ['The data is normalized using the method: "', normFieldName, '"']; i = i + 1;
            headerRow{i} = ['The EEG channel is: "', chName, '"']; i = i + 1;
            headerRow{i} = ['-------']; i = i + 1;
            headerRow{i} = ['The stat field used is: "', statField, '"']; i = i + 1;
            headerRow{i} = ['-------']; i = i + 1;        
            headerRow{i} = ['Subjects done (correspond to the rows below in the matrix)']; i = i + 1;
            for ii = 1 : length(uniqueSubjects)
                uniqueSubjectsString{ii} = [uniqueSubjects{ii}, ', '];
            end        
            headerRow{i} = uniqueSubjectsString; i = i + 1;
            headerRow{i} = ['-------']; i = i + 1;
            headerRow{i} = ['This data is written disk: ', dateStr, ' (year-month-day)']; i = i + 1;
            headerRow{i} = ['===================================================']; i = i + 1;
            
        elseif strcmp(dataType, 'AUX')
            % stimType, chName
            i = 1;
            headerRow{i} = ['Numerical data of the "', dataType, '"']; i = i + 1;
            headerRow{i} = ['===================================================']; i = i + 1;
            headerRow{i} = ['Power computation type is "', stimType, '"']; i = i + 1;            
            headerRow{i} = ['Power band is = "', chName, '"']; i = i + 1;
            headerRow{i} = ['-------']; i = i + 1;
            headerRow{i} = ['Subjects done (correspond to the rows below in the matrix)']; i = i + 1;
            headerRow{i} = ['-------']; i = i + 1;
            headerRow{i} = ['This data is written disk: ', dateStr, ' (year-month-day)']; i = i + 1;
            headerRow{i} = ['===================================================']; i = i + 1;
            
        end

        % define delimiter
        delimiterType = '';

        % write the headers to disk
        for i = 1 : length(headerRow)        
            if i == 1
                dlmwrite(fileOutWithPath, headerRow{i}, 'delimiter', delimiterType)
            else
                dlmwrite(fileOutWithPath, headerRow{i}, '-append', 'delimiter', delimiterType)
            end
        end

    %% DATA
    
        dataColumns = noOfSessions * noOfConditions;
        
        %% data header
        delimiterType = ',';
        j = 1;
        
            % 1st row, the Condition
            for i = 1 : dataColumns
                if i == 1 || i == 5 || i == 9
                    dataHeader{1}{i} = [conditions{j}, delimiterType];
                    j = j + 1;
                else
                    dataHeader{1}{i} = [' ', delimiterType];
                end
            end

            % 2nd row, the Session
            for i = 1 : dataColumns
                
                condit = floor(i/noOfSessions);
                ses = (i - (condit*noOfSessions));
                if ses == 0
                    ses = 4;
                end                  
                dataHeader{2}{i} = [sessions{ses}, delimiterType];
                
                % add the subject header manually
                if i == dataColumns
                    % dataHeader{2}{i+1} = 'SUBJECT';
                    % easier to write numbers later than string
                end
                
            end

            % write to disk
            for rows = 1 : length(dataHeader)
                dlmwrite(fileOutWithPath, dataHeader{rows}, '-append', 'delimiter', '')                
            end
    
    
        %% actual data
        
            % define the data matrix to be written (3D input to 2D output)
        
            for cond = 1 : noOfConditions                
                for ses = 1 : noOfSessions                    
                    ind = (cond - 1) * noOfSessions + ses;                    
                    for sub = 1 : length(uniqueSubjects)                    
                        dataMatrix(sub, ind) = inputMatrix(cond, ses, sub);                        
                    end                    
                end                
            end
    
            dlmwrite(fileOutWithPath, dataMatrix, '-append', 'delimiter', delimiterType)                           
            
            % write averages of the inputMatrix also to the disk (this
            % should be the same as the means stored in statsOut)
            dlmwrite(fileOutWithPath, '  ', '-append', 'delimiter', '')
            dlmwrite(fileOutWithPath, 'Subject averages (+/- SD):', '-append', 'delimiter', '')
            
            averOfSubjects = nanmean(dataMatrix);
            sdOfSubjects = nanstd(dataMatrix);
            dlmwrite(fileOutWithPath, averOfSubjects, '-append', 'delimiter', delimiterType)
            dlmwrite(fileOutWithPath, sdOfSubjects, '-append', 'delimiter', delimiterType)
            
            
            
        %% Write calculated stats?

            % stored there the results
                % stat_results
                % stat_assumptions