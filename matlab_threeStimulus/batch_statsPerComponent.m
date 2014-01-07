function [statsOut, matricesSessionNorm, outlierOut] = ...
    batch_statsPerComponent(dataOut, statsPer, erpBandType, erpComponent, erpFilterType, fieldValue, fileNameFields, fileNames, stimulusType, subjects, outlierFilenameList, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatsOut.mat';
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
    
    %{
    dataOut
    dataOut.bright
    dataOut.bright.session1
    dataOut.bright.session1.target
    dataOut.bright.session1.target.(erpTypes{1})
    dataOut.bright.session1.target.(erpTypes{1}).bl
    dataOut.bright.session1.target.(erpTypes{1}).bl.(erpFilterType)
    dataOut.bright.session1.target.(erpTypes{1}).bl.(erpFilterType).Cz
    %}

    dataOutField = fieldnames(dataOut);   

    % create "intensityBased-matrices (4D)
    % dataNormField
    
    noOfSessions = length(fieldnames(dataOut.(dataOutField{1})));
    erpTypes = fieldnames(dataOut.(dataOutField{1}).session1.target);
    subjectsIn = fieldnames(dataOut.(dataOutField{1}).session1.target.(erpTypes{1}));
    noOfSubjects = length(subjectsIn);    
    noOfTrials = length(dataOut.(dataOutField{1}).session1.target.(erpTypes{1}).(subjectsIn{1}).(erpFilterType).Cz);

    disp('  .. . create intensity matrix (4D)')
    matricesIntensity = batch_getPerSubjectMatrix_Intensity(dataOut, erpBandType, erpComponent, erpFilterType, fieldValue, noOfSubjects, noOfTrials, handles);   

    % Each session have only one value, i.e. the average of 40
    % trials 
    if strcmp(statsPer, 'session')            

        % get only one mean value out from session
        disp('  .. . . average trials together (3D)')
        [matricesSessionAver, outlierOut] = batch_averageSessionTrialsTogether(matricesIntensity, fileNames, subjects, outlierFilenameList, handles);     

        % Normalize using two methods
        matricesSessionNorm.absolute = matricesSessionAver;
        matricesSessionNorm.darkCondition = batch_normalizeComponents(matricesSessionAver, 'darkCondition', subjects, handles);
        matricesSessionNorm.firstSession = batch_normalizeComponents(matricesSessionAver, 'firstSession', subjects, handles);
        % matricesSessionNorm.firstSession = matricesSessionAver;

        disp('  .. . . . calculate stats of these')
        normFieldNames = fieldnames(matricesSessionNorm);
        for i = 1 : length(normFieldNames)

            % simple stats like mean, SD, median
            statsOut.(normFieldNames{i}) = batch_calculateStatsFromStats(matricesSessionNorm.(normFieldNames{i}), handles);

            % calculate significances between different conditions
            % statsOut.(normFieldNames{i}) = batch_calculateStatSignificances(statsOut.(normFieldNames{i}), matricesSessionNorm.(normFieldNames{i}), handles);

        end
        

    % Now a vector is created that has a length of (number of
    % trials per session * number of sessions, e.g. 40 x 4 =
    % 160 trials)
                
    
    elseif strcmp(statsPer, 'trials')
    
        % not necessarily functioning correctly as it has not been
        % used that extensively
        warning(['  !  ', statsPer, ' -based analysis might not work properly, not tested that well!'])

        % concacenate the session together
        disp('  .. . . concatenate sessions together (3D)')
        matricesConcatenated = batch_concatenateSessionsTogether(matricesIntensity, handles);  

        % finally, average the subjects
        disp('  .. . . . calculate stats of these')
        statsOut = batch_averageTheSubjectsMatrix(matricesConcatenated, statsPer, handles);    

    end
    %}
