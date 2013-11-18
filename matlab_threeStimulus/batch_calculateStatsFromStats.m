function sessionStat = batch_calculateStatsFromStats(matricesSessionAver, handles)

    %{
    matricesSessionAver
    matricesSessionAver.dim
    a = matricesSessionAver.dim.Cz    
    b = matricesSessionAver.dim.Cz{1}
    whos
    %}
    ERPtypes = fieldnames(matricesSessionAver);

    noOfSessions = length(matricesSessionAver.(ERPtypes{1}).dim.Cz);
    noOfSubjects = length(matricesSessionAver.(ERPtypes{1}).dim.Cz{1}.mean);

    conditions = fieldnames(matricesSessionAver.(ERPtypes{1}));
    noOfConditions = length(conditions);

    channels = fieldnames(matricesSessionAver.(ERPtypes{1}).(conditions{1}));
    noOfChannels = length(channels);           

    % now we need a column vector with as many rows as there are
    % sessions (i.e. 4) and the subjects are averaged
    dim = 2;
    flag = 0;

    for j = 1 : length(ERPtypes)
        for condition = 1 : noOfConditions
            for ch = 1 : noOfChannels

                % re-calculate the SD                

                for session = 1 : noOfSessions

                    % should be the size of 4 sessions x 2 subjects                    
                    dataMatrix(session, :) = matricesSessionAver.(ERPtypes{j}).(conditions{condition}).(channels{ch}){session}.mean;
                    dataMatrix_SD(session, :) = matricesSessionAver.(ERPtypes{j}).(conditions{condition}).(channels{ch}){session}.SD;

                end

                matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}) = ...
                    batch_calculateStatsForMatrix(dataMatrix, dim, flag, handles);

                % Data in
                % dataMatrix;

                % Mean out
                meanOut = matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}).mean;

                % n of non-NaN values
                nOut = matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}).n;

                % Now we have two SDs:
                % 1) One coming from the averaging of the all the trials together                
                % 2) Another from averaging the subjects                
                SD_1 = dataMatrix_SD;
                    % this one will have rows as many as there are sessions
                    % (i.e. 4) and as many columns as there subjects

                SD_2 = matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}).SD;
                    % only once column as subjects are averaged together,
                    % and as many rows as there sessions

                % What is the best way to get a nice SD estimate:
                SD_1_1col = nanmean(SD_1,2); %?

                    % convert possible NaNs to zeroes
                    SD_2(isnan(SD_2)) = 0;
                    SD_1_1col(isnan(SD_1_1col)) = 0;

                % RMS of the two column vectors?
                SD_rms = sqrt( (SD_1_1col .^ 2) + (SD_2 .^ 2) );
                
                matricesStat.(ERPtypes{j}).(conditions{condition}).(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset}).SD_rms = SD_rms;

                % debug display
                if ch == 2 && 1 == 2
                    disp(conditions{condition})
                    disp(handles.parameters.BioSemi.chName{ch+handles.parameters.BioSemi.chOffset})
                    disp(['Mean  SD_rms   nOut   SD_1_1col SD_2  SD_1_1col(s)'])
                    debug = [meanOut SD_rms nOut SD_1_1col SD_2 SD_1 ];
                    disp(debug)
                end


            end           
        end
    end

    sessionStat = matricesStat;
