function [epochsOut, artifactIndices] = pre_artifactFASTER_wrapper(epochsIn, fixedIndices, EEGfixedIndices, EOGfixedIndices, rawInput, parameters, erpType, handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempFASTERwrapper.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
        else
            if handles.flags.saveDebugMATs == 1
                if ~strcmp(erpType, 'standard')
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
    
    disp(['     (', erpType, ')'])
    
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
                    set(fig, 'Position', [0.10*scrsz(3) 0.075*scrsz(4) 0.82*scrsz(3) 0.85*scrsz(4)])
                    rows = 4;
                    cols = 3;
                    
                sp_i = 1;
                sp(sp_i) = subplot(rows,cols, [1 4 7]);
                    yOffset = 55; % [uV]
                    %plot(1,1)
                    plot_allTheEpochsToSingleSubplot(t*1000, EEG.data(1:parameters.EEG.nrOfChannels+1,:,:), parameters, yOffset)
                    leg(1) = legend('Base', 'Cz', 'Fz', 'Pz', 'Oz', 'EOG');
                        set(leg(1), 'Position',[0.0567851959361391 0.469834826427773 0.0598693759071118 0.120380739081747])
                        legend('boxoff')
                        xlabel('Time [ms]')
                        ylabel('Epochs')            
                        title(['Epochs IN (', erpType, ')'])
                        drawnow
                       
                sp_i = sp_i + 1;
                sp(sp_i) = subplot(rows,cols, [10]);
                
                    % get the average waveform
                    averWaveForm = nanmean(EEG.data(1:parameters.EEG.nrOfChannels+1,:,:),3);               
                    plot(t*1000, averWaveForm)
                    xlabel('Time [ms]')
                    ylabel('Amplitude [\muV]')
                    title(['Average Waveform'])
                    
                    yLimsIns = get(gca, 'YLim');                    

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
                disp(['       ... ', num2str(artifactsRemoved_step2), ' epochs rejected (Step 2) from ', num2str(noOfEpochs1), ' epochs'])
            
                
        %% STEP 3: ICA
                
            % Skip re-referencing as we were referencing for ear electrodes            
            if parameters.artifacts.FASTER_skipICA == 0 || parameters.artifacts.useADJUST == 1
                
                % Do ICA (the runica of EEGLAB, infomax)
                chans_to_interp = 0; % we skipped the step 1            
                k_value = 25;
                ica_chans = 1:4;           
                lpf_band = [w_l-(t_l/2) w_l+(t_l/2)];
                blinkCh = parameters.EEG.nrOfChannels + 1;

                disp('        ... computing ICA (runica), might take some time (try to switch to fastICA for speed)')      
                [EEG, indelec_st3, zs_st3, num_pca, activData, blinkData] = faster_step3_ICA(EEGmatrix, EEG, k_value, ica_chans, chans_to_interp, lpf_band, blinkCh, epochLength, parameters);
                
            else
            
                disp('        ... skipping ICA for FASTER')
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
            epochPerChannelIsArtifacted = zeros(noOfEpochs1,length(eeg_chans));
            epochPerChannelStep2Corrected = zeros(noOfEpochs1,length(eeg_chans));
            
            zs_st4 = zeros(noOfEpochs2, length(eeg_chans), 4); % 4 is the number of different "artifact scores"
            
            for v = 1 : noOfEpochs2
                
                list_properties = single_epoch_channel_properties(EEG, v, eeg_chans);                
                rejection_options.measure = ones(1,size(list_properties,2)); % values of the statistical parameters (see flow chart)
                rejection_options.z = parameters.artifacts.FASTER_zThreshold_step4 * ones(1,size(list_properties,2)); % Z-score threshold
                
                [indelec_st4, zs_st4_raw] = min_z_mod(list_properties,rejection_options);
                zs_st4(v,:,:) = zs_st4_raw;
                indicesOut = logical(min_z(list_properties, rejection_options));
                lengths_ep{v} = eeg_chans(indicesOut);       
                
                % re-correct for baseline
                EEGrecorr = pre_removeBaseline_epochByEpoch(EEG.data(:,:,v)', v, parameters, handles);
                EEG.data(:,:,v) = EEGrecorr';          
                
            end
            
            % no we need to add the artifacted epochs back to the
            % matrix             
            EEGout(:,:,step2_linearIndices_nonArtifact) = EEG.data(:,:,:);                  
            EEG.data = EEGout;
            
            artifactsRemoved_step4 = sum(sum(epochPerChannelIsArtifacted == 1));
            disp(['         ... ', num2str(artifactsRemoved_step4/parameters.EEG.nrOfChannels), ' epochs rejected (Step 4) from ', num2str(noOfEpochs2), ' epochs'])             
            artifactIndices = epochPerChannelIsArtifacted + epochPerChannelStep2Corrected;
            
            % average over channels
            zs_st4 = nanmean(zs_st4, 3);
            
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
                
                disp('          ... using ADJUST for artifact removal')      
                
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
                
                disp('          ... skipping ADJUST for artifact removal')      
                 
            end
            
            
        
        %% CORRECT USING FIXED THRESHOLDS
            
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
            
            artifactsRemoved_fixed = sum(sum(fixedIndices));
            disp(['          ... ', num2str(artifactsRemoved_fixed/parameters.EEG.nrOfChannels), ' epochs rejected (Fixed) from ', num2str(noOfEpochs1), ' epochs'])
        
            % update output
            artifactIndices = logical(artifactIndices + fixedIndices);
            

        %% PLOT ARTIFACT REMOVAL STEPS
        
            if debugFASTER == 1
                if strcmp(erpType, 'target') || strcmp(erpType, 'distracter')  || strcmp(erpType, 'standard')
                    subplotIndices = [2 5 8 11];
                    sp_i = plot_FASTER_steps(fig, sp, sp_i, leg, rows, cols, indelec_st2, zs_st2, indelec_st3, zs_st3, num_pca, activData, blinkData, zs_st4, ...
                        epochPerChannelIsArtifacted, epochPerChannelStep2Corrected, fixedIndices, subplotIndices, parameters, handles);
                end
            end       
        
        
        %% PLOT "ARTIFACT-FREE" EPOCHS and the AVERAGE WAVEFORM
        EEG_mat = pre_epochsVectorToMatrix(EEG.data, epochLength);
        
        if debugFASTER == 1
            if strcmp(erpType, 'target') || strcmp(erpType, 'distracter')  || strcmp(erpType, 'standard')
                                   
                sp_i = sp_i + 1;
                sp(sp_i) = subplot(rows,cols, [3 6 9]);
                    plot_allTheEpochsToSingleSubplot(t*1000, EEG_mat(1:parameters.EEG.nrOfChannels,:,:), parameters, yOffset)
                    leg(2) = legend('Base', 'Cz', 'Fz', 'Pz', 'Oz');
                        set(leg(2), 'Position',[0.911647314949202 0.575867861142218 0.0511611030478955 0.120380739081747])
                        legend('boxoff')
                        xlabel('Time [ms]')
                        ylabel('Epochs')            
                        title(['Epochs OUT (', erpType, ')'])
                        drawnow
                       
                sp_i = sp_i + 1;
                sp(sp_i) = subplot(rows,cols, [12]);
                
                    % get the average waveform
                    averWaveForm = nanmean(EEG_mat(1:parameters.EEG.nrOfChannels,:,:),3);
                    plot(t*1000, averWaveForm)
                    xlabel('Time [ms]')
                    ylabel('Amplitude [\muV]')
                    title(['Average Waveform'])
                    ylim(yLimsIns) % use the same y-limits as for the input to make comparison easier
                    
                set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
                set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
                
            end
            
            %% Auto-SAVE
            try
                if handles.figureOut.ON == 1                     
                    drawnow
                    dateStr = plot_getDateString(); % get current date as string          
                    fileNameOut = sprintf('%s%s%s%s', 'debug_FASTER', '_', erpType, '_', strrep(handles.inputFile, '.bdf', ''), '.png');
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
        
        
%% SUBFUNCTIONS
%% -------------------------------

    function sp_i = plot_FASTER_steps(fig, sp, sp_i, leg, rows, cols, indelec_st2, zs_st2, indelec_st3, zs_st3, num_pca, activData, blinkData, zs_st4, epochPerChannelIsArtifacted, epochPerChannelStep2Corrected, isNaN_fixed, subplotIndices, parameters, handles)

        % STEP 2
        sp_i = sp_i + 1;
        sp(sp_i) = subplot(rows,cols, subplotIndices(1));

            hold on
            plot(zs_st2,'linewidth',1)

            % Thresholds
            plot(1:length(indelec_st2),parameters.artifacts.FASTER_zThreshold*ones(1,length(indelec_st2)),'k-.','linewidth',2)
            plot(1:length(indelec_st2),-parameters.artifacts.FASTER_zThreshold*ones(1,length(indelec_st2)),'k-.','linewidth',2)
            leg(2) = legend('mean dev. from ch.', 'Ep Var', 'Max \DeltaAmplitude', 'Location', 'EastOutside');
            legend('boxoff')
            xlabel('Epochs')
            ylabel(['Z-score (thr = ', num2str(parameters.artifacts.FASTER_zThreshold), ')'])
            title('STEP 2: EPOCHS')
            grid off
            axis tight

        % STEP 3
        %{
        sp_i = sp_i + 1;
        sp(sp_i) = subplot(rows,cols, subplotIndices(2));
        
            plot(activData, blinkData, 'ko', 'MarkerSize', 2)
            xlabel('ICA Activation')
            ylabel('EOG channel ("blinks")')            
            title('STEP 3: Corrcoef()')
            grid off
            axis square
        %}
        
        sp_i = sp_i + 1;
        sp(sp_i) = subplot(rows,cols, subplotIndices(2));       
                 
            if ~isnan(num_pca)
                x = linspace(1,num_pca,num_pca);

                hold on
                if size(zs_st3,3) > 1
                    for ica = 1 : size(zs_st3,1)
                        y = squeeze(zs_st3(ica,:,:));                
                        plot(x,y,'o')
                    end
                else
                    plot(x,zs_st3,'o')
                end

                % Thresholds
                plot(1:length(indelec_st3),parameters.artifacts.FASTER_zThreshold*ones(1,length(indelec_st3)),'k-.','linewidth',2)
                plot(1:length(indelec_st3),-parameters.artifacts.FASTER_zThreshold*ones(1,length(indelec_st3)),'k-.','linewidth',2)
                leg(3) = legend('Median gradient', 'Mean slope', 'Kurtosis', 'Hurst', 'Blink', 'Location', 'EastOutside');
                legend('boxoff')
                xlabel('ICA components')
                ylabel(['Z-score (thr = ', num2str(parameters.artifacts.FASTER_zThreshold), ')'])   
                title('STEP 3: ICA')
                grid off
                axis tight
                set(gca, 'XTick', x, 'XLim', [min(x)-0.5 max(x)+0.5])      
            else
                plot(0, 0)                
                leg(3) = legend('nothing plotted');                    
                    legend('hide')                    
                    axis off
                    title('STEP 3: ICA')
            end
            

        % STEP 4: SINGLE-CHANNEL, SINGLE-EPOCH ARTIFACTS    
        sp_i = sp_i + 1;
        sp(sp_i) = subplot(rows,cols, subplotIndices(3));
        
            hold on
            plot(zs_st4,'linewidth',1)

            % Thresholds
            plot(1:length(indelec_st2),parameters.artifacts.FASTER_zThreshold_step4*ones(1,length(indelec_st2)), 'k-.', 'linewidth', 2)
            plot(1:length(indelec_st2),-parameters.artifacts.FASTER_zThreshold_step4*ones(1,length(indelec_st2)), 'k-.', 'linewidth', 2)
            leg(4) = legend('Variance', 'Median slope', 'Ampl. range', 'Electr. Drift', 'Location', 'EastOutside');
            legend('boxoff')
            xlabel('Epochs')
            ylabel(['Z-score (thr = ', num2str(parameters.artifacts.FASTER_zThreshold_step4), ')'])
            title('STEP 4: SINGLE-Ch, SINGLE-Ep (abs mean of chs)')
            grid off
            % axis tight
            xlim([1 length(zs_st4)])
            
        
        % CONCLUSION
        sp_i = sp_i + 1;
        sp(sp_i) = subplot(rows,cols, subplotIndices(4));
        
            hold on
            b(1) = bar(logical(sum(epochPerChannelStep2Corrected,2)), 'r');
            b(2) = bar(logical(sum(epochPerChannelIsArtifacted,2)), 'g');
            b(3) = bar(logical(sum(isNaN_fixed,2)), 'b');
            hold off
            
            % change alpha of bars (transparency)
            alpha = .75;
            ch = get(b(1),'child');
                set(ch,'facea',alpha)
            ch = get(b(2),'child');
                set(ch,'facea',alpha)
            ch = get(b(3),'child');
                set(ch,'facea',alpha)
            
                xlabel('Epochs')
                ylabel('Artifact (ON/OFF)')            
                title('Artifacted Epochs')
                set(gca, 'XLim', [1 length(isNaN_fixed)], 'YLim', [0 1.2]) 
                
                leg(5) = legend(['Step2, n=', num2str(sum(indelec_st2 == 1), '%3.0f')],...
                                ['Step4, n=', num2str(sum(sum(epochPerChannelIsArtifacted))/parameters.EEG.nrOfChannels, '%3.2f')],...
                                ['Fixed, n=', num2str(sum(sum(isNaN_fixed))/parameters.EEG.nrOfChannels, '%3.2f')]);
                    set(leg(5), 'Position', [0.576378809869376 0.255879059350504 0.0576923076923077 0.0506718924972004])
                    legend('boxoff')
                    

        % STEP 5: GRAND AVERAGE (skip)
        %sp(5) = subplot(rows,cols, subplotIndices(4));

        set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
        set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
       

    function plot_allTheEpochsToSingleSubplot(t, EEGepochs, parameters, yOffset)
       
        numberOfEpochs = size(EEGepochs,3);

        hold on
        for ij = 1 : numberOfEpochs

           yOff = yOffset*ij; % update the horizontal line (baseline)

           y(ij,:,:) = EEGepochs(:,:,ij);

           % Horizontal line (baseline)
           l(ij) = line([min(t) max(t)], [yOff yOff], 'Color', 'k');                                     

           % Filtered
           p_ERP(ij, :) = plot(t, yOff + squeeze(y(ij,:,:)));

           yTickPos(ij) = yOff; % save for yTick positions (Trial)

           yTickLabel{ij} = num2str(ij);
           drawnow

        end

        set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabel)
        yLims = get(gca, 'YLim');        
        yMin = nanmin(nanmin(nanmin(y)));
        yMax = nanmax(nanmax(nanmax(yOff + y(:,:,:))));
        % not optimal for incoming data as gets the min and max from whole
        % matrix, fix later
        try
            set(gca, 'YLim', [yMin yMax])
        catch err            
            if strcmp(err.identifier, 'MATLAB:hg:propswch:PropertyError')
                warning('yLimits are the same, most likely all the epochs are considered as artifacts?')
            end
        end

    function [EEG, indelec_st3, zs_st3, num_pca, activData, blinkData] = faster_step3_ICA(EEGmatrix, EEG, k_value, ica_chans, chans_to_interp, lpf_band, blinkCh, epochLength, parameters)
            
        % quick'n'dirty
        EEG.dataIN = EEG.data;
        EEG.data = EEGmatrix';
        
        num_pca = min(floor(sqrt(size(EEG.data(:,:),2) / k_value)),(size(EEG.data,1) - length(chans_to_interp) - 1));
        num_pca = min(num_pca,length(setdiff(ica_chans,chans_to_interp)));
        
        if num_pca == 0
           warning(['num_pca = ', num2str(num_pca)]) 
           num_pca = 1;
        end

        try
            
            % Original infomax implementation (slow)                
            [EEG.icaweights, EEG.icasphere, compvars, bias, signs, lrates, EEG.icaact] = runica(EEGmatrix', 'extended', 1, 'pca', num_pca, 'verbose', 'off');
            unmixing_matrix = EEG.icaweights*EEG.icasphere;

            % We could use FastICA instead, suggested also in the discussion
            % of FASTER, http://research.ics.aalto.fi/ica/fastica/
            %{
            [A, unmixing_matrix] = fastica(EEGmatrix', 'lastEig', num_pca, 'verbose', 'off', 'displayMode', 'off'); % gives only the estimated mixing matrix A and the separating matrix W.
            
                %size(EEG.icaweights), % number of PCAs x number of channels
                %size(EEG.icasphere), % number of PCAs x number of ch
                %size(unmixing_matrix)
                
                EEG.icaweights = A'; % is this correct?
                %EEG.icasphere ?
                
                % how to define the number of PCAs, and EXTENDED?
                % check that 'lastEig' is the same as above for runica
                

            % compute ICA activation waveforms = weights*sphere*(data-meandata)
            % Usage: >> [activations] = icaact(data,weights,datamean);
            %}
            
            % EEGLAB variables, see e.g. http://sccn.ucsd.edu/wiki/A05:_Data_Structures                
                EEG.icachansind = ica_chans;
                EEG.trials = length(EEGmatrix) / epochLength;
                EEG.pnts = epochLength;                
                
            EEG.srate = parameters.EEG.srate;
            EEG.icaact = icaact(EEGmatrix', unmixing_matrix);
            EEG.icawinv = pinv(EEG.icaweights*EEG.icasphere); % http://sccn.ucsd.edu/pipermail/eeglablist/2009/002907.html                

            % size(EEG.icaact) % number of PCAs x dataSamples
            
            
        catch err
            err
            error('improve error catching!!')
        end

        % after the ICA routine we have the EEG data as 2-dimensional
        % matrix, and we need to to separate the epochs to the third
        % dimension
        [list_properties, activData, blinkData] = component_properties_mod(EEG, blinkCh,lpf_band);
        rejection_options.measure=ones(1,size(list_properties,2)); % values of the statistical parameters (see flow chart)
        rejection_options.z = parameters.artifacts.FASTER_zThreshold * ones(1,size(list_properties,2)); % Z-score threshold
        [indelec_st3, zs_st3] = min_z_mod(list_properties,rejection_options); % rejected components
        
        step3_linearIndices = find(indelec_st3);
        if ~isempty(step3_linearIndices)
            disp(['          - subtracting ICA artifacts, found ', num2str(length(step3_linearIndices)), ' artifacted ICA activation channels']) 
            for i = 1 : length(step3_linearIndices)
                for ch = 1 : ica_chans
                    EEG.dataIN(:,ch) = EEG.icaact(:, step3_linearIndices(i));
                end
            end
        else
            disp(['          - No ICA artifacts found'])             
        end
        
        % quick'n'dirty
        EEG.data = EEG.dataIN;
        