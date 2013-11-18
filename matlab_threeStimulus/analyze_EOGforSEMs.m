function SEM = analyze_EOGforSEMs(EOG, EOG_raw, handles);

    % EOG  - bandpass filtered with "general filter"
    %        [parameters.filter.bandPass_loFreq parameters.filter.bandPass_hiFreq]
    % EOG  - as recorded originally with no manipulations, just trimmed to
    %        match the trigger ON for recording trigger
    
    % handles  - all the settings, parameters, etc.
    
    SEM.dummy = [];