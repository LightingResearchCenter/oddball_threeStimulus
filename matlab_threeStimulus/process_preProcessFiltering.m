function [dataMatrix_filtGeneral, firstBandpassFilteredMatrix, artifactNaN_indices, dataMatrix_filtAlpha, dataMatrix_filt, dataMatrix_filt_CNV] = process_preProcessFiltering(dataMatrixIn, handles)

    % [rowsIn, colsIn] = size(dataMatrixIn); % e.g. a lot of rows x 8 channels
    offset = 2; % get rid of reference channels

    % remove the baseline (reference)
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
        [dataMatrix_filtGeneral, firstBandpassFilteredMatrix, artifactNaN_indices] ...
            = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_loFreq, handles.parameters.filter.bandPass_hiFreq, handles.parameters.filterOrder, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings                 


     %% ALPHA                   

        % Alpha-filtered ERP, as done by Barry et al. (2000) for
        % example, http://dx.doi.org/10.1016/S0167-8760(00)00114-8                   
        dataType = 'Alpha';
        matrixIn = firstBandpassFilteredMatrix;
        artifactIndices = artifactNaN_indices;
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filtAlpha,~,~] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_Alpha_loFreq, handles.parameters.filter.bandPass_Alpha_hiFreq, handles.parameters.filterOrder_Alpha, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings


    %% ERP

        % For ERP (P3, N2, P3-N2)
        dataType = 'ERP';
        matrixIn = firstBandpassFilteredMatrix;
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filt,~,~] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_ERP_loFreq, handles.parameters.filter.bandPass_ERP_hiFreq, handles.parameters.filterOrder, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings


    %% CNV

        % For CNV
        dataType = 'CNV';
        matrixIn = firstBandpassFilteredMatrix;
        [rowsIn, colsIn] = size(matrixIn);
        rawMax = max(max(abs(matrixIn(:,1:handles.parameters.EEG.nrOfChannels))));

        % use subfunction wrapper
        [dataMatrix_filt_CNV,~,~] = pre_componentArtifactFiltering(matrixIn, artifactIndices, rawMax, rowsIn, colsIn, ... % EEG matrix parameters
            handles.parameters.filter.bandPass_CNV_loFreq, handles.parameters.filter.bandPass_CNV_hiFreq, handles.parameters.filterOrder_CNV, handles.parameters.filterOrderSteep, ... % filter parameters
            dataType, handles.parameters, handles); % general settings


        % Debug plot for bandpass characteristics
        if handles.flags.showDebugPlots == 1
            plot_bandPassFilter(dataMatrix_filt_CNV(:,1:handles.parameters.EEG.nrOfChannels),...
                                dataMatrix_filt(:,1:handles.parameters.EEG.nrOfChannels),...
                                dataMatrix(:,1+offset:handles.parameters.EEG.nrOfChannels+offset), handles)
        end

