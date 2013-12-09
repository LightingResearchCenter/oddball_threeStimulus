function dataOut = pre_artifactRemovalInputData(dataIn, EOG, ECG, i, dataType, parameters, handles)

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempArtifactRemoval.mat';
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
    %}
        
    [dataSamples, dataChannels] = size(dataIn);
    dataOut = dataIn;
    
    %% Use Muscle detection with inverse filtering (BioSig)
    
        % do only once, with the general filtering
        if strcmp(dataType, 'General') == 1 
            
            % https://github.com/donnchadh/biosig/blob/master/biosig/t250_ArtifactPreProcessingQualityControl/detectmuscle.m
            
            % References: 
            % [1]	Van de Velde, M., Van Erp, G., Cluitmans, P., 1998. 
            %	Detection of muscle artefact in the normal human awake EEG. 
            % 	Electroencephalography and Clinical Neurophysiology 107 (2), 149-158.
            %   http://dx.doi.org/10.1016/S0013-4694(98)00052-2

            [rows,cols] = size(dataIn);

            iter = 1;
            Mode = 1;

            S = zeros(rows,cols);          
            E_muscle = zeros(rows,cols);
            
            for i = 1 : cols
                [~, E_muscle(:,i)] = detectmuscle(dataIn(:,i), iter, Mode);
                %[INI{i},S(:,i), E_muscle(:,i)] = detectmuscle(dataIn(:,i), iter, Mode);
            end

            %{
            whos        
            hold on
            plot(S(:,1),'r')
            plot(E(:,1), 'g')
            plot(dataIn(:,1), 'b')
            hold off        
            pause
            %}

            numberOfNaNs = length(dataOut(isnan(E_muscle) == 1));
            S(~isnan(E_muscle),:) = NaN;
            NaNPercentage = (numberOfNaNs / (rows*cols)) * 100;
            disp(['     . Muscle artifact detection algorithm - ',  'Number of NaNs: ', num2str(numberOfNaNs), ', percentage: ', num2str(NaNPercentage), '%'])           

            dataOut_afterMuscle = dataIn;
            dataOut_afterMuscle(isnan(E_muscle)) = NaN;
            dataOut = dataOut_afterMuscle;
            
        else
            
            disp(['     . No muscle artifacts removed (only for "general" condition)'])
            dataOut_afterMuscle = dataIn;
            numberOfNaNs = 0;
            
        end


    %% EOG Correction

        % Requires BioSig
        % http://sourceforge.net/projects/biosig/

        % REGRESS_EOG yields the regression coefficients for 
        % correcting EOG artifacts in EEG recordings.
        % http://biosig.sourceforge.net/help/biosig/t250/regress_eog.html

        % combine EOG and EEG
        S1 = [dataOut EOG];

        % find the artifacts
        EL = [1 2 3 4];
        OL = 5;
        try
            R = regress_eog(S1, EL, OL);
        catch err
            if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                error('EOG Rejection not done! Have you installed the BioSig toolbox needed for EOG artifact rejection? http://sourceforge.net/projects/biosig/')
            else
                err
            end                
        end

        % apply the correction, with offset correction
        S2 = [ones(size(S1,1),1),S1] * R.r1;
        % S2b = S1 * R.r0;    % without offset correction
        dataOut = S2(:,EL);
        
        % save value for debug
        dataOut_afterEOG_removal = dataOut;       
            

    %% ECG Rejection

        % - You can expect the best performance when R = regress_eog(..) is 
        % computed from a data segment with large EOG and EKG artifacts (ask the 
        % subject to perform different eye movements at the beginning of the 
        % recording).
    
        % Same logic as for EOG rejection above, just use the ECG
        % channel as the OL
        % combine EOG and EEG
        S1 = [dataOut ECG];

        % find the artifacts
        EL = [1 2 3 4];
        OL = 5;
        try
            R = regress_eog(S1, EL, OL);
        catch err
            if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                error('ECG Rejection not done! Have you installed the BioSig toolbox needed for EOG artifact rejection? http://sourceforge.net/projects/biosig/')
            else
                err
            end  
            
        end

        % apply the correction, with offset correction
        S2 = [ones(size(S1,1),1),S1] * R.r1;
        % S2b = S1 * R.r0;    % without offset correction
        dataOut = S2(:,EL);

        % save value for debug
        dataOut_afterECG_removal = dataOut;       

    %% Alternatively one could use ICA

        % not used now: [Eeglablist] Using ICA with EEG, ECG and EOG.
        % http://sccn.ucsd.edu/pipermail/eeglablist/2012/004965.html
        % ---> Schlögl A, Ziehe A, Müller K-R. 2009. Automated ocular artifact removal: comparing regression and component-based methods. Nature Precedings. http://dx.doi.org/10.1038/npre.2009.3446.1.
        % http://precedings.nature.com/documents/3446/version/1/files/npre20093446-1.pdf        
        
        % we leave the option to do this after the epoch-by-epoch correction
        % see PROCESS_singleFile the general pipeline
        
        
    %% DEBUG
    
        if handles.flags.showDebugPlots == 1
            
            close all
            
            t = linspace(1, length(dataIn(:,1)), length(dataIn(:,1))) / parameters.EEG.srate;
            ch = 2; % 1 for Fz, 2 for Pz  
           
            zoomRange = [160 165]; % in seconds
            plot_artifactRemovalPlot(t, dataIn, dataOut_afterMuscle, dataIn, dataOut_afterEOG_removal, ...
                                        dataOut_afterECG_removal, dataOut, ch, zoomRange, handles.parameters, handles)
            
        end


