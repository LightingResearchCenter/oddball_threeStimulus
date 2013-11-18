function dataOut_ICA = pre_artifactByICA(dataIn, i, parameters, handles)

    debugMatFileName = 'tempICA.mat';
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

    % EEGLAB Tutorial
    % Data Analysis
    % I.9. Performing Independent Component Analysis of EEG data
    % http://cognitrn.psych.indiana.edu/busey/temp/eeglabtutorial4.301/maintut/ICA_decomposition.html    
    
    dataOut_ICA = dataIn;
    regular = dataIn.oddball_regular;
    irregular = dataIn.oddball_irregular;
    sampleVector = (linspace(1,length(regular),length(regular)))';       
    
    %% PRE-PROCESS
    
        % Now the data cannot contain any NaN nor Inf values, so remove them
        isnanVector_regular = pre_removeNaNs(regular);
        isnanVector_irregular = pre_removeNaNs(irregular);

            % trim the NaNs out now
            sampleVector_reg_woNaNs = sampleVector(~isnanVector_regular);
            regular_woNaNs = regular(~isnanVector_regular,:);
            sampleVector_irreg_woNaNs = sampleVector(~isnanVector_irregular);
            irregular_woNaNs = irregular(~isnanVector_irregular,:);

            % display    
            disp(['        .. regular epochs had ', num2str(100*(1-(length(regular_woNaNs)/length(isnanVector_regular)))), '% of NaNs'])
            disp(['        .. irregular epochs had ', num2str(100*(1-(length(irregular_woNaNs)/length(isnanVector_irregular)))), '% of NaNs'])

        % if you want the ICA progress to be displayed on command window (hogs
        % up the command window rather nicely)
        if parameters.artifacts.show_ICA_verbose == 1
            verboseFlag = 'on';
        else
            verboseFlag = 'off';
        end

    %% RUNICA()
        
        % Apply the runica() function
        try 
            [weights.regular,sphere.regular] = runica(regular_woNaNs', 'stop', 1E-7, 'verbose', verboseFlag);
            [weights.irregular,sphere.irregular] = runica(irregular_woNaNs', 'stop', 1E-7, 'verbose', verboseFlag);        
        catch err        
            if strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
               error('The data contains either NaN or Inf which is not allowed for the SVD done by ICA algorithm') 
            else
                err
            end        
        end        
    
    %% Now you could remove data based on the ICA 
    
        % see e.g. http://www.nbtwiki.net/doku.php?id=tutorial:automatic_and_semi-automatic_methods_for_eeg_pre-processing#.UfFvIt9jEx0
        % there are EEGLAB PLUGINs for automatic rejection such as:
        % * FASTER (Fully Automated Statistical Thresholding for EEG artifact Rejection)
        %   http://dx.doi.org/10.1016/j.jneumeth.2010.07.015
        %   http://www.mee.tcd.ie/neuraleng/Research/Faster
        
        % * ADJUST (Automatic EEG artifact detection based on the joint use of spatial and temporal features)
        %   http://dx.doi.org/10.1111/j.1469-8986.2010.01061.x
        %   http://www.unicog.org/pm/pmwiki.php/MEG/RemovingArtifactsWithADJUST
        
        warning('   NOTHING ACTUALLY DONE currently with ICA decomposition. Implement this!')
        
            % FASTER looked more promising
    
    %% PLOT
    
        % plot_artifactRemovalPlot_ICA(dataIn, EOG, ECG, dataOut, dataOut_ICA, parameters, handles) 
        whos
        pause
    
    
    function isnanVector = pre_removeNaNs(matrixIn)        
        
        isnanMatrix = isnan(matrixIn);

        % Now add the indices together so that the vector sizes match after NaN
        % removal if the NaN indices are not the same for each vector
        [rows,cols] = size(isnanMatrix);    
        isnanVector = zeros(rows,1);
        for i = 1 : cols
            isnanVector = isnanVector + isnanMatrix(:,i);
        end

        % convert back to logical
        isnanVector = logical(isnanVector);
    
    
    
    function plot_artifactRemovalPlot_ICA(dataIn, EOG, ECG, dataOut, dataOut_ICA, parameters, handles)
        
        scrsz = handles.style.scrsz;
        fig = figure('Color', 'w');
            set(fig, 'Position', [0.35*scrsz(3) 0.08*scrsz(4) 0.55*scrsz(3) 0.85*scrsz(4)])
            
            
    
       