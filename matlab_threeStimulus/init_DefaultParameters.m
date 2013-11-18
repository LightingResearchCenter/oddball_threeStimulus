function parameters = init_DefaultParameters(handles)

    %% BIOSEMI Channel Definitions
    
        parameters.BioSemi.chName{1} = 'Ref_RightEar'; 
        parameters.BioSemi.chName{2} = 'Ref_LeftEar'; 
        parameters.BioSemi.chName{3} = 'Fz'; % ch3 - EX3: Fz
        parameters.BioSemi.chName{4} = 'Cz'; % ch4 - EX4: Cz
        parameters.BioSemi.chName{5} = 'Pz'; % ch5 - EX5: Pz
        parameters.BioSemi.chName{6} = 'Oz'; % ch6 - EX6: Oz
        parameters.BioSemi.chName{7} = 'EOG'; % ch7 - EX7: EOG (was put below the right eye)
        parameters.BioSemi.chName{8} = 'HR'; % ch8 - EX8: Heart Rate (was put on the chest)   
        parameters.BioSemi.chOffset = 2; % number of channels to omit (i.e. the reference channels)

            % The input trigger signals are saved in an extra channel (the status channel), with the same sample rate as the electrode channels.
            % http://www.biosemi.com/faq/trigger_signals.htm
            parameters.BioSemi.chName{9} = 'Status/Trigger';            

        % Inverse polarity (if P300 is negative and N2 is positive)
        parameters.EEG.invertPolarity = 0;

        % Defines the signals to be read from the BDF file   
        parameters.SignalsToRead = [1 2 3 4 5 6 7 8 9]; % Channels 1-8 and 9 is the Trigger

        % Define triggers

            % Button 1 & 2
            % Standard Tone
            % Oddball Tone
            % Irregular Cycle
            % Recording start 

            % http://www.biosemi.com/faq/trigger_signals.htm
            parameters.triggerSignals.buttons = [1 2];
            parameters.triggerSignals.oddTone = 9;
            parameters.triggerSignals.stdTone = 10;            
            parameters.triggerSignals.distracterTone = 11;
            parameters.triggerSignals.recON = 12;  
            
            parameters.triggerSignals.audioON = 5;

            parameters.triggerPrecision = 24; % 24 bits, http://www.biosemi.com/faq/trigger_signals.htm

        % General EEG Parameters
        % parameters.EEG.srate - read automatically during import of individual BDF file
        parameters.EEG.nrOfChannels = 4; % number of EEG channels
            
    %% Band-pass filter parameters
    
        % for discussion of the filter limits, you can see for example:
        
            % Acunzo DJ, MacKenzie G, van Rossum MCW. 2012. 
            % Systematic biases in early ERP and ERF components as a result of high-pass filtering. 
            % Journal of Neuroscience Methods 209:212–218. 
            % http://dx.doi.org/10.1016/j.jneumeth.2012.06.011.
            
            % Widmann A, Schröger E. 2012. 
            % Filter effects and filter artifacts in the analysis of electrophysiological data. 
            % Frontiers in Perception Science:233. 
            % http://dx.doi.org/10.3389/fpsyg.2012.00233.
    
        % GENERAL
        parameters.filter.bandPass_loFreq = 0.01;
        parameters.filter.bandPass_hiFreq = 50;
        parameters.filterOrder = 6; % filter order   
        parameters.filterOrderSteep = 100; % additional pass of filtering with steeper cut
        parameters.applySteepBandPass = 0;
        
            % 0.01 Hz recommended as low-cut off frequency for ERP by:
            % * Acunzo et al. (2012), http://dx.doi.org/10.1016/j.jneumeth.2012.06.011
            % * Luck SJ. 2005. An introduction to the event-related potential technique. Cambridge, Mass.: MIT Press.
        
        % ALPHA, see Barry et al. (2000), http://dx.doi.org/10.1016/S0167-8760(00)00114-8        
        parameters.filter.bandPass_Alpha_loFreq = 8;
        parameters.filter.bandPass_Alpha_hiFreq = 13;
        parameters.filterOrder_Alpha = 10; % filter order   
    
        % ERP
        parameters.filter.bandPass_ERP_loFreq = parameters.filter.bandPass_loFreq;
        parameters.filter.bandPass_ERP_hiFreq = 20;
        parameters.filterOrder_ERP = parameters.filterOrder; % filter order   
    
        % parameters for re-bandbass filtering for extracting the CNV
        parameters.filter.bandPass_CNV_loFreq = parameters.filter.bandPass_loFreq;
        parameters.filter.bandPass_CNV_hiFreq = 6;
        parameters.filterOrder_CNV = parameters.filterOrder; % filter order   
        
    
    %% Artifact rejection parameters    
    
        % Fixed thresholds
        parameters.artifacts.fixedThr = 150; % fixed threshold (uV) of artifacts    
                                             % 100 uV in Molnar et al. (2008), http://dx.doi.org/10.1111/j.1469-8986.2008.00648.x
        parameters.artifacts.fixedThrEOG = 150; % fixed threshold (uV) of EOG artifacts
                                               % 70 uV in e.g. Acunzo et al. (2012), http://dx.doi.org/10.1016/j.jneumeth.2012.06.011
        parameters.artifacts.applyFixedThrRemoval = 1; % convert values above threshold to NaN
        parameters.artifacts.applyFixedThrEOGRemoval = 1; % convert values above threshold to NaN
                        
            % fixed threshold "detrending parameters"
            % not that crucial if some distortion is introduced by too
            % aggressive low cut of 1 Hz for example, but without
            % detrending, finding threshold exceeding samples would be
            % rather impossible due to possible DC drifts and trends
            parameters.artifacts.fixedDetrendingLowCut = parameters.filter.bandPass_loFreq;
            parameters.artifacts.fixedDetrendingHighCut = 300;
            parameters.artifacts.fixedDetrendingOrder = parameters.filterOrder;
        
        % "Advanced artifact removal"
        parameters.artifacts.applyRegressEOG = 1; % apply regress_eog from BioSig to eliminate EOG/ECG-based artifacts
        parameters.artifacts.epochByEpochRemoveBaseline = 0; % use rmbase() to remove baseline before ICA
        parameters.artifacts.useICA = 0; % otherwise use the regress_eog
        parameters.artifacts.show_ICA_verbose = 1;
    
    %% Power spectrum analysis parameters
    
        parameters.powerAnalysis.tukeyWindowR = 0.10; % 0.10 equals to 10% cosine window
        parameters.powerAnalysis.segmentLength = 4.0; % xx second segment lengths for PWELCH
        parameters.powerAnalysis.nOverlap = 50; % overlap [%] between successive segments
        
        chOffset = 2; %EX1 and EX2 for earlobes
        parameters.powerAnalysis.alphaRange = [parameters.filter.bandPass_Alpha_loFreq parameters.filter.bandPass_Alpha_hiFreq];
        parameters.powerAnalysis.alphaCh = [5 6] - chOffset; % Pz and Oz (-2 for ref channels)
                
        parameters.powerAnalysis.eegBins.freqs{1} = [-1.5 1.5];
        parameters.powerAnalysis.eegBins.label{1} = 'alpha'; % on Pz and Oz
        parameters.powerAnalysis.eegBins.ch{1} = [5 6] - chOffset; % from what channel(s) are calculated

        parameters.powerAnalysis.eegBins.freqs{2} = [-1.5 0];
        parameters.powerAnalysis.eegBins.label{2} = 'lowAlpha'; % on Pz and Oz
        parameters.powerAnalysis.eegBins.ch{2} = [5 6] - chOffset; 
        
        parameters.powerAnalysis.eegBins.freqs{3} = [0 1.5];
        parameters.powerAnalysis.eegBins.label{3} = 'highAlpha'; % on Pz and Oz
        parameters.powerAnalysis.eegBins.ch{3} = [5 6] - chOffset; 

        parameters.powerAnalysis.eegBins.freqs{4} = [13 30];
        parameters.powerAnalysis.eegBins.label{4} = 'beta'; % Fz and Cz
        parameters.powerAnalysis.eegBins.ch{4} = [3 4] - chOffset; 

        parameters.powerAnalysis.eegBins.freqs{5} = [5 7];
        parameters.powerAnalysis.eegBins.label{5} = 'theta'; % Fz and Cz
        parameters.powerAnalysis.eegBins.ch{5} = [3 4] - chOffset; 
        
        % Bands used in Lockley et al. 2006 (check abstract)
        % http://www.ncbi.nlm.nih.gov/pubmed/16494083
        parameters.powerAnalysis.eegBins.freqs{6} = [0.5 5.5];
        parameters.powerAnalysis.eegBins.label{6} = 'deltaThetaLockley';
        parameters.powerAnalysis.eegBins.ch{6} = [3 4] - chOffset;
        
        % Total band power, used for calculating the ratios
        parameters.powerAnalysis.eegBins.freqs{7} = [0.5 parameters.filter.bandPass_hiFreq];
        parameters.powerAnalysis.eegBins.label{7} = 'totalFzCz';
        parameters.powerAnalysis.eegBins.ch{7} = [3 4] - chOffset;
        
        parameters.powerAnalysis.eegBins.freqs{8} = [0.5 parameters.filter.bandPass_hiFreq];
        parameters.powerAnalysis.eegBins.label{8} = 'totalOzPz';
        parameters.powerAnalysis.eegBins.ch{8} = [5 6] - chOffset;
        
            % for more details see: http://www.mathworks.com/help/signal/ref/pwelch.html
            

    %% EEG FRACTAL ANALYSIS
    
    
    %% ECG ANALYSIS       
    
        parameters.heart.downsampleFactor = 8;
    
        parameters.heart.detrendMethod = 'ougp'; % 'linear' / 'constant' /'ougp'
        parameters.heart.filterWindow = 'hanning'; % 'hanning' / 'hamming' / 'blackman' / 'bartlett'
        parameters.heart.noOfFFTPoints = 'auto'; % 'auto' / '512' / '1024'
        parameters.heart.spectralAlgorithm = 'FFT'; % 'FFT' / 'AR model'
        parameters.heart.rrTimeSeriesUpsampleRate = 4; % [Hz], see e.g. http://dx.doi.org/10.1007/s11517-012-0928-2

        parameters.heart.freqBins.ULF = [0.0005 0.003]; % [Hz]
        parameters.heart.freqBins.ULFStar = [0.002 0.01]; % [Hz]
        parameters.heart.freqBins.VLF = [0.003 0.01]; % [Hz]
        parameters.heart.freqBins.LF = [0.04 0.15]; % [Hz]
        parameters.heart.freqBins.HF = [0.15 0.40]; % [Hz]
        parameters.heart.freqBins.TOTAL = [0.005 2]; % [Hz]

        parameters.heart.freqResolution = 0.005;
        parameters.heart.freqRange = [parameters.heart.freqBins.TOTAL(1) parameters.heart.freqBins.TOTAL(2)+parameters.heart.freqResolution];
        parameters.heart.ougpBinByBin = 0;

        if parameters.heart.freqRange(1) == 0
            warning('Low cutoff for HRV analysis can''t be zero, use the freqResolution for example! Script will work but PSD will be full of NaNs')
        end
        
            
    %% Time-Frequency Analysis
    
        parameters.timeFreq.logFrq = 0; % whether you use linear or log frequency scale
                
        % Scales is the parameter a
        % scale factor a is defined as the inverse of frequency (a=1/ƒ)
        % in general. ERPWaveLab wants the scales to be a list of
        % frequencies that you want to analyze       
        parameters.timeFreq.scaleLimits = [1 30];
        parameters.timeFreq.numberOfFreqs = 2*(parameters.timeFreq.scaleLimits(2) - parameters.timeFreq.scaleLimits(1)) + 1;
        parameters.timeFreq.tukeyWindowR = 0.05;
                
        parameters.timeFreq.bandwidthParameter = 2;
        parameters.timeFreq.waveletResol = 12;
        parameters.timeFreq.timeResolutionDivider = 16;
            % e.g. with 8,  freqRange: [3.8721 495.6228], freqResolution: 0.0422
            % e.g. with 16, freqRange: [1.9360 247.8114], freqResolution: 0.0211
            % e.g. with 32, freqRange: [0.9680 123.9057], freqResolution: 0.0105. timeRes: 7.8125 ms
        

            
    %% 3-Stimulus Oddball parameters
    
        parameters.oddballTask.numberOfTrialsPerCycle = 8;
        parameters.oddballTask.repeatsOfCycle = 40; % e.g. 8 x 40 = 320 trials
        parameters.oddballTask.nrOfStandardsPerCycle = 6; % 75% of the tones in other words
        
        parameters.oddballTask.targetFreq = 1000;
        parameters.oddballTask.distractFreq = [500 5000];
        parameters.oddballTask.standardFreq = 500;
        
        parameters.oddballTask.triggerDuration = 0.2; % [s]
        parameters.oddballTask.SOA_duration = 2.0; % [s], 
        parameters.oddballTask.ERP_duration = 0.5; % 0.50; % [s]
        parameters.oddballTask.ERP_baseline = 0.5; % [s], needed e.g. for ep_den
            % in seconds, note that the trigger is for the whole 800 ms, from
            % which we can trim some of the end away
            
            % DC Offset / Detrend correction (pre_epochToERPs --> )
            %{
            parameters.epochERP.detrendingLowCut = parameters.artifacts.fixedDetrendingLowCut;
            parameters.epochERP.detrendingHighCut = parameters.artifacts.fixedDetrendingHighCut;
            parameters.epochERP.detrendingOrder = 4;
            %}
            
            % Epoch baseline removel (remove the mean of the pre-stimulus
            % baseline period for example)
            parameters.oddballTask.baselineRemove_index1 = 1;
            parameters.oddballTask.baselineRemove_index2 = 0 - (-parameters.oddballTask.ERP_baseline); % 500 ms

            % for pre-stimulus power analysis
            % see e.g. Barry et al. (2000), http://dx.doi.org/10.1016/S0167-8760(00)00114-8
            parameters.oddballTask.preERP_power_segmentLength = parameters.oddballTask.ERP_baseline; % [s]
            parameters.oddballTask.preERP_power_tukeyWindowR = parameters.powerAnalysis.tukeyWindowR;
            parameters.oddballTask.preERP_power_nOverlap = parameters.powerAnalysis.nOverlap;
            parameters.oddballTask.preERP_IAF_range = [-4 1.5]; % from Klimesch 1999, or more tight [-3.5 1.0]

        % Fixed time windows for the ERP components (see. Jongsma 2006,
        % 2013), http://dx.doi.org/10.1016/j.clinph.2006.05.012 and
        % http://dx.doi.org/10.1016/j.clinph.2012.09.009

            % -300 to 0 ms both in Jongsma 2006 and 2013
            parameters.oddballTask.timeWindows.CNV = [-0.3 0];
                    % [-110 -10] ms in Molnar et al. (2008), http://dx.doi.org/10.1111/j.1469-8986.2008.00648.x

            % 150 ms to 250 ms in Min 2013 (http://dx.doi.org/10.1016/j.brainres.2013.09.031)
            parameters.oddballTask.timeWindows.N1 = [0.150 0.250];            
                                    
            % 180 ms to 220 ms in Jongsma 2006, 180ms to 280 ms in Jongsma 2013
            parameters.oddballTask.timeWindows.N2 = [0.180 0.270]; 
                % 400 to 600 ms in Min 2013

            % 350 ms to 430 ms in Jongsma 2006, 280 to 480 ms in Jongsma 2013
            parameters.oddballTask.timeWindows.P3 = [0.280 0.480];

    %% EP_DEN Parameters
    
        % parameters.ep_den.sr = 512; % sampling rate
        % parameters.ep_den.stim = 513; % stim
        % parameters.ep_den.samples = 1024; % number of samples

        % Jongsma et al. denoised the epoch twice, first for extracting the CNV
        % with a different scale setting, and then second time for the
        % remaining components    
        parameters.ep_den.scales_postStim = 8; % number of scales
        parameters.ep_den.scales_preStim = 10; % number of scales

        parameters.ep_den.plot_type ='coeff';  % 'coeff', 'bands', 'single', 'contour'
        parameters.ep_den.den_type = 'do_den'; %' do_den' or 'load_den_coeff' 
        parameters.ep_den.auto_den_type = 'NZT';  % 'Neigh' or 'NZT'

        % Sigmoid fit parameters
        parameters.sigmoid.sigmoidFunc = ''; % 'param4'

    %% PLOTTING
    
        handles.parameters.plot.ComponentAveragedSamples = 10; % for component trend plot
        handles.parameters.plot.timeFreq_contourLevels = 64;
        

    %% BATCH ANALYSIS
    
        % number of different intensities
        handles.parameters.batch.noOfIntensities = 3;
        handles.parameters.batch.noOfSessions = 4;
