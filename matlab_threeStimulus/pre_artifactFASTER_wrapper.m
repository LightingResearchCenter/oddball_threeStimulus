function [epochsOut, artifactIndices] = pre_artifactFASTER_wrapper(epochsIn, fixedIndices, EEGfixedIndices, EOGfixedIndices, ...
                                            NaN_indices_moving, NaN_indices_movingEOG, NaN_indices_step, ...
                                            referenceInput, vDiffOutMovWindow, vDiffOutMovWindowEOG, vDiffOutStep, ...
                                            parameters, erpType, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempFASTERwrapper.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
        else
            if handles.flags.saveDebugMATs == 1
                if ~strcmp(erpType, 'fff') % 'standard')
                    % do not save for standard tone as there are so many
                    % trials that debugging and developing of this function
                    % is so much slower compared to target and distracter
                    path = handles.path;
                    save('debugPath.mat', 'path')
                    save(fullfile(path.debugMATs, debugMatFileName))            
                end
            end
        end 
    end
    
    disp(['        - ', erpType])    
    epochsOut = epochsIn;    
    debugFASTER = 1;    
    
    %EEGfixedIndices
    %EOGfixedIndices
    
    parameters.EEG.downsampleFactor = parameters.EEG.srate / parameters.EEG.downsampledRate;
    parameters.EEG.srate = parameters.EEG.downsampledRate;
    
    w_l = parameters.artifacts.FASTER_lpf_wFreq;
    t_l = parameters.artifacts.FASTER_lpf_bandwidth;

    parameters.BioSemi.channelLocationsFile = 'centralChannelLocations.ced'; % could open when importing file
    parameters.artifacts.FASTER_skipICA = 1;
    parameters.artifacts.useADJUST = 0;
    
    %% For a quick introduction to FASTER, see:
        
        % http://www.nbtwiki.net/doku.php?id=tutorial:automatic_and_semi-automatic_methods_for_eeg_pre-processing#.Uouqah8wnNA

        % For more detailed info, see

        % Anthony J Bell and Terrence J Sejnowski (1995). An information-maximisation approach to blind separation and blind deconvolution. 
        % Neural Computation. November 1995, Vol. 7, No. 6, Pages 1129-1159.
         
        % Nolan H, Whelan R, Reilly RB.(2010), FASTER: Fully Automated Statistical Thresholding for EEG artifact Rejection. 
        % J Neurosci Methods. 2010 Sep 30;192(1):152-62. Epub 2010 Jul 21.
        
        % MATLAB CODE: http://www.mee.tcd.ie/neuraleng/Research/Faster

    %% Actual workflow
    
        % get the matrix of the epochs
        EEGmatrix_orig = epochsIn.ERP(:, 1:parameters.EEG.nrOfChannels+2);
        
            % epochsIn (e.g. for standard tone)
            %                RT: [1x240 double]
            %   samplesPerEpoch: 4096
            %           Indices: {1x240 cell}
            %               ERP: [983040x6 double]
            % whos            
        
            inputType = 'matrix';
            epochLength = epochsIn.samplesPerEpoch;
        
            % downsample
            t_orig = linspace(-parameters.oddballTask.ERP_baseline, parameters.oddballTask.ERP_duration, epochLength);
            t = linspace(-parameters.oddballTask.ERP_baseline, parameters.oddballTask.ERP_duration, epochLength/parameters.EEG.downsampleFactor);
            
            x = linspace(1,length(EEGmatrix_orig), length(EEGmatrix_orig));
            x_i = linspace(1,length(EEGmatrix_orig), length(EEGmatrix_orig) / parameters.EEG.downsampleFactor);

            epochLength = epochLength / parameters.EEG.downsampleFactor;
            epochsOut.samplesPerEpoch = epochLength; % to return structure
            
            % Downsample            
            for ch = 1 : size(EEGmatrix_orig,2)
                EEGmatrix(:,ch) = interp1(x, EEGmatrix_orig(:,ch), x_i);
            end
            
        % now the epochs are concatenated to a single vector but the EEGLAB
        % version is to have 3D matrix where the 3rd dimension is the number of
        % epochs
        EEG_mat = pre_epochsVectorToMatrix(EEGmatrix',epochLength);
        EEG.data = EEG_mat;


        %% PLOT EPOCHS
        if debugFASTER == 1            

            if strcmp(erpType, 'target') || strcmp(erpType, 'distracter')  || strcmp(erpType, 'standard')

                scrsz = [1 1 1680 1050];
                fig = figure('Color','w','Name',[erpType, ': FASTER Debug']);
                    set(fig, 'Position', [0.02*scrsz(3) 0.075*scrsz(4) 0.95*scrsz(3) 0.85*scrsz(4)])
                    rows = 4;
                    cols = 4;
                    
                noOfValidTrials = sum(~isnan(permute(EEG.data(1:parameters.EEG.nrOfChannels,:,:),[2 1 3])),3); % per channel
                noOfValidTrials = noOfValidTrials(1,:);
                    
                sp_i = 1;
                sp(sp_i) = subplot(rows,cols, [1 5 9]);
                    yOffset = 55; % [uV]
                    %plot(1,1)
                    try
                        plot_allTheEpochsToSingleSubplot(sp(sp_i),t*1000, EEG.data(1:parameters.EEG.nrOfChannels+1,:,:), parameters, yOffset)                    
                        leg(1) = legend('Base', 'Cz', 'Fz', 'Pz', 'Oz', 'EOG');
                            set(leg(1), 'Position',[0.0567851959361391 0.469834826427773 0.0598693759071118 0.120380739081747])
                            legend('boxoff')
                            lab(1,1) = xlabel('Time [ms]');
                            lab(1,2) = ylabel('Epochs');
                            tit(1) = title(['Epochs IN (', erpType, ')']);
                            xlim([min(t*1000) max(t*1000)])
                            drawnow
                    catch err
                        err
                        % not robust enough if you are looking for example
                        % the TARGET plot while Matlab is plotting for
                        % standard, then you might get an invalid handle
                        % error
                    end
                       
                sp_i = sp_i + 1;
                sp(sp_i) = subplot(rows,cols, [13]);
                
                    % get the average waveform
                    averWaveForm = nanmean(EEG.data(1:parameters.EEG.nrOfChannels+1,:,:),3);       
                    hold on
                    line([min(t*1000) max(t*1000)], [0 0], 'Color', 'k')
                    plot(t*1000, averWaveForm)
                    hold off                    
                    lab(2,1) = xlabel('Time [ms]');
                    lab(2,2) = ylabel('Amplitude [\muV]');
                    tit(2) = title(['Average Waveform, n = [', num2str(noOfValidTrials), ']']);
                    xlim([min(t*1000) max(t*1000)])
                    yLimsIns = get(gca, 'YLim');                    
                    drawnow

            end
        end

        %% STEP 1: CHANNELS (skip)

            % Skip this for our data

            % For the explanation of steps, see Fig. 2 of the original
            % article, same there: 
            % http://www.nbtwiki.net/lib/exe/detail.php?id=tutorial%3Aautomatic_and_semi-automatic_methods_for_eeg_pre-processing&media=tutorial:faster_flow_chart.jpg

        %% STEP 2: EPOCHS

            % this list_properties matrix has as many rows as there are epochs, and
            % three columns corresponding to different parameters that FASTER uses
            % to determine whether the epoch has an artifact or not:
            % column 1: Epoch's mean deviation from channel means
            % column 2: Epoch variance
            % column 3: Max amplitude difference
            eeg_chans = 1:parameters.EEG.nrOfChannels;
            list_properties = epoch_properties_mod(EEG, EEGmatrix', eeg_chans, inputType, epochLength);

            rejection_options.measure = ones(1,size(list_properties,2)); % values of the statistical parameters (see flow chart)
            rejection_options.z = parameters.artifacts.FASTER_zThreshold * ones(1,size(list_properties,2)); % Z-score threshold

            %a = rejection_options.measure
            %b = rejection_options.z

            % rejected epochs (returns indices and z-values)
            [indelec_st2, zs_st2] = min_z_mod(list_properties,rejection_options);
            noOfEpochs1 = length(indelec_st2);
            
            % indelect_st2 is a logical vector for artifacted epochs, we
            % can get linear indices to aid the combining of the artifact
            % indices with the STEP 4
            step2_linearIndices = find(indelec_st2);
            step2_linearIndices_nonArtifact = ~indelec_st2;
            
            % Actually reject (exclude the NaN
            EEG.data = EEG.data(:,:,indelec_st2 == 0);
                artifactsRemoved_step2 = sum(indelec_st2 == 1);
                disp(['         ... ', num2str(artifactsRemoved_step2), ' epochs rejected (Step 2) from ', num2str(noOfEpochs1), ' epochs'])
            
                
        %% STEP 3: ICA
                
            % Skip re-referencing as we were referencing for ear electrodes            
            if parameters.artifacts.FASTER_skipICA == 0 || parameters.artifacts.useADJUST == 1
                
                % Do ICA (the runica of EEGLAB, infomax)
                chans_to_interp = 0; % we skipped the step 1            
                k_value = 25;
                ica_chans = 1:4;           
                lpf_band = [w_l-(t_l/2) w_l+(t_l/2)];
                blinkCh = parameters.EEG.nrOfChannels + 1;

                disp('          ... computing ICA (runica), might take some time (try to switch to fastICA for speed)')      
                [EEG, indelec_st3, zs_st3, num_pca, activData, blinkData] = pre_FASTER_step3_ICA(EEGmatrix, EEG, k_value, ica_chans, chans_to_interp, lpf_band, blinkCh, epochLength, parameters);
                
            else
            
                disp('          ... skipping ICA for FASTER')
                indelec_st3 = NaN;
                zs_st3 = NaN;
                num_pca = NaN;
                activData = NaN;
                blinkData = NaN;
                
                
            end

        %% STEP 4: CHANNELS IN EPOCHS
        
            noOfEpochs2 = size(EEG.data,3);
            lengths_ep=cell(1,noOfEpochs2);
            
            EEGout = zeros(size(EEG.data,1), size(EEG.data,2), noOfEpochs1);
            EEGout(:,:,:) = NaN;
            %epochPerChannelIsArtifacted = zeros(noOfEpochs1,length(eeg_chans));
            %epochPerChannelStep2Corrected = zeros(noOfEpochs1,length(eeg_chans));
            
            zs_st4_allCh = zeros(noOfEpochs2, length(eeg_chans), 4); % 4 is the number of different "artifact scores"
            
            for v = 1 : noOfEpochs2
                
                %{
                subplot(2,1,1)                
                plot(EEG.data(:,:,v)')
                %}
                
                list_properties = single_epoch_channel_properties(EEG, v, eeg_chans);                
                rejection_options.measure = ones(1,size(list_properties,2)); % values of the statistical parameters (see flow chart)
                rejection_options.z = parameters.artifacts.FASTER_zThreshold_step4 * ones(1,size(list_properties,2)); % Z-score threshold
                
                [indelec_st4, zs_st4_raw] = min_z_mod(list_properties,rejection_options);
                    % [noOfChannels, noOfMeasures] = size(zs_st4_raw)
                    % noOfChannels = length(indelec_st4)
                    
                zs_st4_allCh(v,:,:) = zs_st4_raw;
                indicesOut(v,:) = logical(min_z(list_properties, rejection_options));
                lengths_ep{v} = eeg_chans(indicesOut(v,:));       
                
                % actually reject the epoch (convert the artifacted epochs
                % to NaNs)
                linearIndices = indicesOut(v,:) == 1;
                EEG.data(linearIndices,:,v) = NaN;
                
                % re-correct for baseline (useful actually only if you had
                % done ICA subtraction in STEP 3)
                EEGrecorr = pre_removeBaseline_epochByEpoch(EEG.data(:,:,v)', v, parameters, handles);
                EEG.data(:,:,v) = EEGrecorr';          
                
                % combine the 4 different measures so that we have an absolute max
                % per channel, and if this exceeds the z-threshold, then
                % the epoch is artifacted
                [absMeasuresForChannelsPerEpoch,I] = max(abs(zs_st4_raw), [], 2);
                for ch = 1 : length(eeg_chans) % keep the sign, for loop maybe more easier to read by human later                                        
                    measuresForChannelsPerEpoch(ch) = zs_st4_raw(ch,I(ch));
                end
                zs_st4(v,:) = measuresForChannelsPerEpoch;
                %{
                subplot(2,1,2)
                plot(EEG.data(:,:,v)')
                title(['indicesOut: ', num2str(indicesOut(v,:))]);
                legend('Fz', 'Cz', 'Pz', 'Oz', 'EOG', 'ECG')
                pause
                %}
                              
            end
            
            % no we need to add the artifacted epochs back to the
            % original matrix  (as Step 4 is not using the epochs rejected
            % in Step 2)
            EEGout(:,:,step2_linearIndices_nonArtifact) = EEG.data(:,:,:);                  
            EEG.data = EEGout;
            
                % the same for zs_st4
                zs_st4temp = zeros(noOfEpochs1, length(eeg_chans));
                zs_st4temp(:,:) = NaN;
                zs_st4temp(step2_linearIndices_nonArtifact,:) = zs_st4;
                zs_st4 = zs_st4temp;
                
            
            %artifactsRemoved_step4 = sum(sum(epochPerChannelIsArtifacted == 1));
            epochPerChannelIsArtifacted_step4 = logical((squeeze(sum(isnan(EEG.data),2)))');
            epochPerChannelIsArtifacted_step4 = epochPerChannelIsArtifacted_step4(:,1:length(eeg_chans));
                  
                % now this step4 contains step 2 also, separate if you want
            
            artifactsRemoved_step4 = sum(indicesOut == 1);
            disp(['           ... ', num2str(artifactsRemoved_step4), ' epochs rejected (Step 4) from ', num2str(noOfEpochs2), ' epochs'])                       
            
            artifactIndices = epochPerChannelIsArtifacted_step4 + repmat(indelec_st2,1,length(eeg_chans));
            
            
            %% INTERPOLATION (skip)
            
                %a = EEG.srate
                %EEG

                % Interpolate if you want, also interpolates using channel
                % locations for not having discontinuities in the scalp
                % topography which we don't need, we could just convert the
                % epochs found to NaNs

                %{
                % Add EEGLAB fields to avoid crashing
                EEG.setname = erpType;
                EEG.nbchan = length(eeg_chans);
                EEG.xmax = parameters.oddballTask.ERP_duration;             
                EEG.xmin = -parameters.oddballTask.ERP_baseline;
                EEG.chanlocs = [];
                try 
                    EEG = h_epoch_interp_spl_mod(EEG, lengths_ep, blinkCh, handles);
                catch err
                    h = 1
                    err                
                end
                %EEG
                %}
                

        %% STEP 5: GRAND AVERAGE (skip)

        
        %% ADJUST (EEGLAB Plugin / Algorithm) for remaining epochs
        
            if parameters.artifacts.useADJUST == 1
                
                disp('            ... using ADJUST for artifact removal')      
                
                % get the channel locations
                EEG.chanlocs = readlocs(parameters.BioSemi.channelLocationsFile);
                
                % slightly cumbersome now as in our EEG.data the 2 last
                % channels are EOG and ECG
                EEGorig = EEG.data;
                EEG.data = EEG.data(1:parameters.EEG.nrOfChannels,:,:);
                EEG.filename = 'dkk';
                
                % need o come up with the .chanlocs for our chanlocs, in
                % order to make this work
                [art, horiz, vert, blink, disc,...
                    soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
                    soglia_GDSF, GDSF, soglia_V, nuovaV, soglia_D, maxdin] = ADJUST_mod(EEG, 'outStringFileNameDummy', handles)

                % re-assign
                EEG.data = EEGorig;
            
            else                
                disp('            ... skipping ADJUST for artifact removal')
            end            
            
        
        %% REJECT USING FIXED THRESHOLDS
            
            % obtained in "PROCESS_singleFile.m" using the subfunction
            % "pre_artifactFASTER_fixedThresholds.m"
            % fixedIndices          
                       
            for ch =  1 : size(fixedIndices,2)                
                for ep =  1 : size(fixedIndices,1)            
                    if fixedIndices(ep,ch) == 1
                        EEG.data(ch,:,ep) = NaN;
                    end
                end
            end           
            
            %artifactsRemoved_fixed = sum(sum(fixedIndices));            
            artifactsRemoved_fixed = sum(fixedIndices);
            artifactsRemoved_fixedEEG = sum(EEGfixedIndices);            
            artifactsRemoved_fixedEOG = sum(EOGfixedIndices);     
            
            % replicate to all the channels
            if size(NaN_indices_movingEOG,2) == 1
                NaN_indices_movingEOG = repmat(NaN_indices_movingEOG, 1, size(NaN_indices_moving,2));
            end
            
            try                
                artifacts_CRAP = NaN_indices_moving + NaN_indices_step + NaN_indices_movingEOG;                
            catch err              
                size(NaN_indices_moving)
                size(NaN_indices_step)
                size(NaN_indices_movingEOG)
                err.identifier
                error('Dimensions do not match')
                whos
            end
            artifactsRemoved_CRAP = sum(artifacts_CRAP);
            disp(['             ... ', num2str(artifactsRemoved_CRAP), ' epochs rejected (CRAP) from ', num2str(noOfEpochs1), ' epochs'])
            disp(['              ... ', num2str(artifactsRemoved_fixedEEG), ' epochs rejected (Fixed EEG) from ', num2str(noOfEpochs1), ' epochs'])
            disp(['               ... ', num2str(artifactsRemoved_fixedEOG), ' epochs rejected (Fixed EOG) from ', num2str(noOfEpochs1), ' epochs'])
            disp(['                ... ', num2str(artifactsRemoved_fixed), ' epochs rejected (Fixed: EEG+EOG+CRAP) from ', num2str(noOfEpochs1), ' epochs'])
          
            % update output
            %faster_artifactIndices = logical(artifactIndices + fixedIndices);            
            fixedThresholdIndices = EEGfixedIndices + EOGfixedIndices;
            artifactIndices = artifactIndices + fixedThresholdIndices;
            
        %% REJECT USING the CRAP 
        
            %NaN_indices_moving, NaN_indices_step
            for ch =  1 : size(NaN_indices_moving,2)                
                for ep =  1 : size(NaN_indices_moving,1)            
                    if NaN_indices_moving(ep,ch) == 1
                        EEG.data(ch,:,ep) = NaN;
                    end
                    if NaN_indices_movingEOG(ep,ch) == 1
                        EEG.data(ch,:,ep) = NaN;
                    end
                    if NaN_indices_step(ep,ch) == 1
                        EEG.data(ch,:,ep) = NaN;
                    end
                end
            end  
            
            artifactIndices = artifactIndices + NaN_indices_movingEOG + NaN_indices_moving + NaN_indices_step;
            artifactIndices = logical(artifactIndices);

        %% PLOT ARTIFACT REMOVAL STEPS
        
            if debugFASTER == 1
                
                if strcmp(erpType, 'target') || strcmp(erpType, 'distracter')  || strcmp(erpType, 'standard')
                
                    % CRAP and FIXED step
                    subplotIndices_CRAP = [2 6 10 14];                             
                    
                    % condition a bit, only one value per epoch
                    vDiffOutStep = squeeze(nanmax(permute(vDiffOutStep, [3 2 1])));
                    vDiffOutMovWindow = (squeeze(nanmax(permute(vDiffOutMovWindow, [3 2 1]))))';                    
                    vDiffOutMovWindowEOG = (squeeze(nanmax(permute(vDiffOutMovWindow, [3 2 1]))))';     
                    
                    try
                        [sp_i, sp] = plot_CRAPandFIXED_steps(fig, sp, sp_i, leg, rows, cols, ...
                            NaN_indices_moving, NaN_indices_movingEOG, NaN_indices_step, EEGfixedIndices, EOGfixedIndices, vDiffOutMovWindow, vDiffOutMovWindowEOG, vDiffOutStep, ...
                            subplotIndices_CRAP, parameters, handles);                   
                    catch err
                        err
                    end
                        
                    % FASTER steps
                    subplotIndices_FASTER = [3 7 11 15];                    
                    [sp_i, sp] = plot_FASTER_steps(fig, sp, sp_i, leg, rows, cols, ...
                        zs_st2, indelec_st3, zs_st3, num_pca, activData, blinkData, zs_st4, ...
                        epochPerChannelIsArtifacted_step4, indelec_st2, artifacts_CRAP, subplotIndices_FASTER, parameters, handles);                    
                    
                end
            end       
        
        
        %% PLOT "ARTIFACT-FREE" EPOCHS and the AVERAGE WAVEFORM
        EEG_mat = pre_epochsVectorToMatrix(EEG.data, epochLength);
        
        if debugFASTER == 1
            if strcmp(erpType, 'target') || strcmp(erpType, 'distracter')  || strcmp(erpType, 'standard')
                                   
                sp_i = sp_i + 1;
                sp(sp_i) = subplot(rows,cols, [4 8 12]);
                    plot_allTheEpochsToSingleSubplot(sp(sp_i),t*1000, EEG_mat(1:parameters.EEG.nrOfChannels,:,:), parameters, yOffset)
                    leg(2) = legend('Base', 'Cz', 'Fz', 'Pz', 'Oz');
                        set(leg(2), 'Position',[0.911647314949202 0.575867861142218 0.0511611030478955 0.120380739081747])
                        legend('boxoff')
                        lab(3,1) = xlabel('Time [ms]');
                        lab(3,2) = ylabel('Epochs');
                        tit(3) = title(['Epochs OUT (', erpType, ')']);
                        xlim([min(t*1000) max(t*1000)])
                        drawnow
                       
                sp_i = sp_i + 1;
                sp(sp_i) = subplot(rows,cols, [16]);                
                    
                    noOfValidTrials = sum(~isnan(permute(EEG_mat(1:parameters.EEG.nrOfChannels,:,:),[2 1 3])),3); % per channel
                    noOfValidTrials = noOfValidTrials(1,:);
                
                    % get the average waveform
                    averWaveForm = nanmean(EEG_mat(1:parameters.EEG.nrOfChannels,:,:),3);
                    hold on
                    line([min(t*1000) max(t*1000)], [0 0], 'Color', 'k')
                    plot(t*1000, averWaveForm)
                    hold off
                    lab(4,1) = xlabel('Time [ms]');
                    lab(4,2) = ylabel('Amplitude [\muV]');
                    tit(4) = title(['Average Waveform, n = [', num2str(noOfValidTrials), ']']);
                    ylim(yLimsIns) % use the same y-limits as for the input to make comparison easier
                    xlim([min(t*1000) max(t*1000)])
                    
                set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
                set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
                set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
                set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold') 
                
            end
            
            %% Auto-SAVE
            try
                if handles.figureOut.ON == 1                     
                    drawnow
                    dateStr = plot_getDateString(); % get current date as string          
                    fileNameOut = sprintf('%s%s%s%s', 'debug_FASTER', '_', strrep(handles.inputFile, '.bdf', ''), erpType,  '_', '.png');
                    export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                    %cd(path.code)
                end
            catch err
                err
                str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
                error(str)
            end
            
        end
        
        %% RETURN
        
            % now the EEG.data is in the format of:
            % channels (without ECG) x samples per epoch x epoch, and we need to
            % concatenate it back to -> channels x samples in all the
            % epochs total
            EEG_concat = reshape(EEG.data, size(EEG.data, 1), (size(EEG.data, 2) * size(EEG.data, 3)));
            epochsOut.ERP = EEG_concat';

            nrOfEpochsOut = size(EEG.data, 3);
            %epochsOut
            %epochsIn

