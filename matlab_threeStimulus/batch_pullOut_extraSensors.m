function [heart, fractalEEG, eye] = batch_pullOut_extraSensors(fileNameFields, handles)

    %% DEBUG
    debugMatFileName = 'tempDataOutExtra.mat';
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

    numberOfFilesFound = length(fileNameFields);
    filesPerSubject = 12;
    numberOfSubjectsDone = numberOfFilesFound/filesPerSubject;
    disp([num2str(numberOfFilesFound), ' files found, thus most likely ', num2str(numberOfSubjectsDone), ' subjects done?'])
    
    if rem(numberOfSubjectsDone, 1) % if non-integer
       error('Lazy programmer have not managed yet the situation where you do not all the sessions done for each subject')       
    end
    
    % erpComponent
    
    for i = 1 : numberOfFilesFound                       
        
        subject = fileNameFields{i}.subject;
        
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
        disp(['    ... load the data from file: "', fileNameFields{i}.fileName, '"'])
        dataIn = load(fullfile(handles.path.matFilesOut, fileNameFields{i}.fileName));
        
        
        %% HEART
        scalarFields = fieldnames(dataIn.analyzed_extraSensors.heart.scalar);
        vectorFields = fieldnames(dataIn.analyzed_extraSensors.heart.vector);
           
            % scalars
            for field = 1 : length(scalarFields)
                heart.(intensity).(session).(subject).scalar.(scalarFields{field}) = dataIn.analyzed_extraSensors.heart.scalar.((scalarFields{field}));
            end
            
            % vectors
            for field = 1 : length(vectorFields)
                % vector fields take quite a lot of memory and are not
                % really needed at this point, correct this if you want
                % some vectors
                % heart.(intensity).(session).(subject).vector.(vectorFields{field}) = dataIn.analyzed_extraSensors.heart.vector.((vectorFields{field}));
                heart.(intensity).(session).(subject).vector.(vectorFields{field}) = []; 
            end
            
        %% FRACTAL        
        for ch = 1 : handles.parameters.EEG.nrOfChannels
            
            scalarFields = fieldnames(dataIn.analyzed_fractal{ch}.scalar);
            % vectorFields = fieldnames(dataIn.analyzed_fractal{ch}.vector)
            vectorFields = {'dummy'};
            
            % scalars
            for field = 1 : length(scalarFields)
                fractalEEG.(intensity).(session).(subject).scalar.(scalarFields{field}){ch} = dataIn.analyzed_fractal{ch}.scalar.((scalarFields{field}));
            end
            
            % vectors            
            for field = 1 : length(vectorFields)
                % fractalEEG.(intensity).(session).(subject).vector.(vectorFields{field}) = dataIn.analyzed_fractal{ch}.vector.((vectorFields{field}));
                fractalEEG.(intensity).(session).(subject).vector.(vectorFields{field}){ch} = [];
            end
            
        end
            
        %% EYE Movements
        
            %scalarFields = fieldnames(dataIn.analyzed_extraSensors.SEM.scalar)
            scalarFields = {'dummy'};
            %vectorFields = fieldnames(dataIn.analyzed_extraSensors.SEM.vector)
            vectorFields = {'dummy'};
            
            % scalars
            for field = 1 : length(scalarFields)
                eye.(intensity).(session).(subject).scalar.(scalarFields{field}) = [];
                %eye.(intensity).(session).(subject).scalar.(scalarFields{field}) = dataIn.analyzed_extraSensors.SEM.scalar.(scalarFields{field});
            end
            
            % vectors            
            for field = 1 : length(vectorFields)
                eye.(intensity).(session).(subject).vector.(vectorFields{field}) = [];
                %eye.(intensity).(session).(subject).vector.(vectorFields{field}) = dataIn.analyzed_extraSensors.SEM.vector.(vectorFields{field});
            end        
        
        
    end