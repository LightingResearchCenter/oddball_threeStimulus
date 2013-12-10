function [dataMatrix_filtGeneral, dataMatrix_filtGeneralRegress, onlyNotchFilteredMatrix, artifactNaN_indices, dataMatrix_filtAlpha, dataMatrix_filt, dataMatrix_filt_CNV] = process_preProcessFiltering(dataMatrixIn, handles)

    % [rowsIn, colsIn] = size(dataMatrixIn); % e.g. a lot of rows x 8 channels
    offset = 2; % get rid of reference channels

    % remove the baseline (reference)
    disp(['   remove reference (mean of ear electrodes)'])
    dataMatrix = pre_removeReference(dataMatrixIn(:,1), dataMatrixIn(:,2), dataMatrixIn, offset, handles.parameters, handles);
    
    %% GENERAL

        % General filtering if you need or want to compute something
        % general about the data which might be
        % * Power Analysis
        % * Detrended Fluctuation Analysis (DFA)
        % * Fractal Length analysis
        %   e.g. 
        dataType = 'General';
        matrixIn = dataMatrix;
        artifactIndices = []; % get them only once
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filtGeneral, dataMatrix_filtGeneralRegress, onlyNotchFilteredMatrix, artifactNaN_indices] ...
            = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_loFreq, handles.parameters.filter.bandPass_hiFreq, handles.parameters.filterOrder, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings                 


     %% ALPHA                   

        % Alpha-filtered ERP, as done by Barry et al. (2000) for
        % example, http://dx.doi.org/10.1016/S0167-8760(00)00114-8                   
        dataType = 'Alpha';
        matrixIn = onlyNotchFilteredMatrix;
        artifactIndices = artifactNaN_indices;
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filtAlpha,~,~,~] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_Alpha_loFreq, handles.parameters.filter.bandPass_Alpha_hiFreq, handles.parameters.filterOrder_Alpha, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings

    %% P300
    
        % You could test a P300 specific filtering either between 0-4 Hz or
        % 0.001 - 4 Hz (to remove DC bias), this band is however mostly used in 
        % BCI application and this filtering will deform the "typical shape" of P300
        
        % see short discussion in:        
            % Ghaderi F, Kim SK, Kirchner EA. 2014. 
            % Effects of eye artifact removal methods on single trial P300 detection, a comparative study. 
            % Journal of Neuroscience Methods 221:41–47. 
            % http://dx.doi.org/10.1016/j.jneumeth.2013.08.025   
        

    %% ERP

        % For ERP (P3, N2, P3-N2)
        dataType = 'ERP';
        matrixIn = onlyNotchFilteredMatrix;
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filt,~,~,~] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_ERP_loFreq, handles.parameters.filter.bandPass_ERP_hiFreq, handles.parameters.filterOrder, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings


    %% CNV

        % For CNV
        dataType = 'CNV';
        matrixIn = onlyNotchFilteredMatrix;
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filt_CNV,~,~,~] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_CNV_loFreq, handles.parameters.filter.bandPass_CNV_hiFreq, handles.parameters.filterOrder_CNV, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings


        % Debug plot for bandpass characteristics
        if handles.flags.showDebugPlots == 1
            plot_bandPassFilter(dataMatrix_filt_CNV(:,1:handles.parameters.EEG.nrOfChannels),...
                                dataMatrix_filt(:,1:handles.parameters.EEG.nrOfChannels),...
                                dataMatrix(:,1+offset:handles.parameters.EEG.nrOfChannels+offset), handles)
        end

