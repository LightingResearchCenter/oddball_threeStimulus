function EEG = import_sortEEGchannels(dataMatrix)

    % BIOSEMI Channel Definitions
    %{
    parameters.BioSemi.chName{1} = 'Ref_RightEar'; 
    parameters.BioSemi.chName{2} = 'Ref_LeftEar'; 
    parameters.BioSemi.chName{3} = 'Fz'; % ch3 - EX3: Fz
    parameters.BioSemi.chName{4} = 'Cz'; % ch4 - EX4: Cz
    parameters.BioSemi.chName{5} = 'Pz'; % ch5 - EX5: Pz
    parameters.BioSemi.chName{6} = 'Oz'; % ch6 - EX6: Oz
    parameters.BioSemi.chName{7} = 'EOG'; % ch7 - EX7: EOG (was put below the right eye)
    parameters.BioSemi.chName{8} = 'HR'; % ch8 - EX8: Heart Rate (was put on the chest)   
    %}

    EEG.reference = (dataMatrix(:,1) + dataMatrix(:,2)) / 2;    
    EEG.Fz = dataMatrix(:,3);
    EEG.Cz = dataMatrix(:,4);
    EEG.Pz = dataMatrix(:,5);
    EEG.Oz = dataMatrix(:,6);
    EEG.EOG = dataMatrix(:,7);
    EEG.HR = dataMatrix(:,8);