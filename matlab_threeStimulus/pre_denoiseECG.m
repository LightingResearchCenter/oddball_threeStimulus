function ECG = pre_denoiseECG(t, ECG_raw, Fs, parameters, handles)
        
    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempECG_Denoising.mat';
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
    
    % debug param
    timeOffset = 0.8;
    zoomWindow = round([mean(t)-timeOffset mean(t)+timeOffset]);
    debugPlotON = 0;
    
    % Re-assign variables
    ECG = ECG_raw;
    style = handles.style;
    scrsz = handles.style.scrsz;
    
    parameters.heart.skip_emdDenoising = 1;
    parameters.heart.EMD_denoising.M1 = 3;
    parameters.heart.EMD_denoising.IM2 = 2;
    parameters.heart.EMD_denoising.T_mult = 0.7;
    parameters.heart.EMD_denoising.threstype = 'soft';
    parameters.heart.EMD_denoising.nofsifts = 8;
    parameters.heart.EMD_denoising.iterations = 15;
    parameters.heart.EMD_denoising.method = 'CIIT';
    parameters.heart.EMD_denoising.altermethod = 'perm';
        
    parameters.heart.skip_waveletReconstruction = 0;
    parameters.heart.reconstr_Par = 5; % 7 in Kabir and Shahnaz (2012)
    parameters.heart.reconstr_motherWavelet = 'Symmlet';
    
    parameters.heart.skip_waveletSoftThresholding = 1;
        % open example.m from WavDen to see the other option, now just
        % using the 'S'   
    

    %{
    % the EMD algorithm takes superlong time at 4,096 Hz, so we need to downsample the signal   
    parameters.heart.downsampleFactor = = 1
    Fs = Fs / parameters.heart.downsampleFactor;
    x_new = (linspace(1, length(ECG), length(ECG) / parameters.heart.downsampleFactor))';
    ECG = pre_nyquistInterpolation(ECG, length(x_new), 1, 0); % Nyquist interpolation  
    t = linspace(min(t), max(t), length(x_new))';
    %}
      
    % truncate to the smaller next power of 2 (alternatively you could
    % interpolate)
    l = length(ECG);
    nn = nextpow2(l);
    indEnd = 2^(nn-1);
    ECG = ECG(1:indEnd);
    t = t(1:indEnd);
    
    % mean of 90-100% percentile
    %ECG_highest_10percentPercentile = nanmean(prctile(ECG,[95 100]))
    
    %% BANDPASS FILTER
        % remove trends with a bandpass filter, and some high-frequency noise
        % with a range of 2.3-30 Hz, see for example Data acquisition lab notes
        % from MIT Media Lab: http://alumni.media.mit.edu/~emunguia/pdf/arithmia%20detector.pdf
        % and order 200 for fir1 filter
        cutOffs = [30 2.3];
        filterOrder = 10; % conservative 
        try 
            ECG = pre_bandbassFilter(ECG, Fs, cutOffs, filterOrder, filterOrder, handles);
        catch
            err
            filterOrder = filterOrder / 2;
            ECG = pre_bandbassFilter(ECG, Fs, cutOffs, filterOrder, filterOrder, handles);
            error('Bandpass filter parameters too "extreme"')
        end
                
    
    
    
    %% Plot input
    if debugPlotON == 1
        fig = figure('Color','w');
            set(fig, 'Position', [0.05*scrsz(3) 0.05*scrsz(4) 0.90*scrsz(3) 0.80*scrsz(4)])
        rows = 5; cols = 3; 
        indices = [1 2 3]; colorPlot = [0 0.6 1]; 
        titleStr = ['INPUT, f_s = ', num2str(Fs), ' Hz'];
        pre_plotQuickECG(t, ECG, zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
        drawnow
    end
    
    %% DENOISE with EMD DENOISING
        
        parameters.heart.skip_emdDenoising = 1;
        if parameters.heart.skip_emdDenoising ~= 1
            
            % COMPUTATIONALLY VERY INEFFICIENT, TAKES AROUND 33 SEC for the
            % one session ECG downsampled to 128 Hz!
            disp('           + EMD Denoising the ECG Signal (very slow)')          
                
                % NOTE! The EMD parameters are not really optimized
                % currently to denoise ECG, for inspiration and guidance
                % you might want to check the following papers:
                
                % Kabir MA, Shahnaz C. 2012. 
                % Denoising of ECG signals based on noise reduction algorithms in EMD and wavelet domains. 
                % Biomedical Signal Processing and Control 7:481–489. 
                % http://dx.doi.org/10.1016/j.bspc.2011.11.003.
                
                % Li H, Wang X, Chen L, Li E. 2013. 
                % Denoising and R-Peak Detection of Electrocardiogram Signal Based on EMD and Improved Approximate Envelope. 
                % Circuits Syst Signal Process:1–16. 
                % http://dx.doi.org/10.1007/s00034-013-9691-3.
                
            % use the Matlab code provided by Yannis Kopsinis
            % http://www.see.ed.ac.uk/~ykopsini/software.html
            % -------------------------------
            % Kopsinis Y, McLaughlin S. 2009. 
            % Development of EMD-Based Denoising Methods Inspired by Wavelet Thresholding. 
            % IEEE Transactions on Signal Processing 57:1351–1362. 
            % http://dx.doi.org/10.1109/TSP.2009.2013885.        

            % Parameters
            M1 = parameters.heart.EMD_denoising.M1;
            IM2 = parameters.heart.EMD_denoising.IM2;
            T_mult = parameters.heart.EMD_denoising.T_mult;
            threstype = parameters.heart.EMD_denoising.threstype;
            nofsifts = parameters.heart.EMD_denoising.nofsifts;
            iterations = parameters.heart.EMD_denoising.iterations;
            method = parameters.heart.EMD_denoising.method;
            altermethod = parameters.heart.EMD_denoising.altermethod; 
            
            heartrateSampleRate = parameters.EEG.srate/parameters.heart.downsampleFactor;

            % Call the subfunction
            whos
            tic
            [ECG_denoised, IMF] = EMDdenoise_modPT(ECG,method,iterations,altermethod,nofsifts,threstype,T_mult,M1,IM2);
            emd_denoisingTime = toc;                

            if debugPlotON == 1           
                indices = [4 5 6]; colorPlot = [0.8 0.2 1]; 
                titleStr = ['EMD Denoising, t = ', num2str(emd_denoisingTime, '%3.1f')];
                pre_plotQuickECG(t, ECG_denoised, zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
                drawnow
            end
            
        else
            disp('           - Skipping the EMD Denoising for the ECG Signal')
        end

    %% DENOISE WITH WAVELET RECONSTRUCTION
    
        if parameters.heart.skip_waveletReconstruction ~= 1
            
            disp(['           + Wavelet Reconstruction of ECG Signal (', parameters.heart.reconstr_motherWavelet, num2str(parameters.heart.reconstr_Par), ')'])
            % Now we can wavelet denoise the ECG and we have a problem with
            % finding the best mother wavelet, Rahman and Riheen suggest the 
            % "Biorthogonal (bior3.9)" to be the best for ECG denoising and
            % "sym5" best from Symmlet family, see:

                % Rahman MW, Riheen MA. 2013. 
                % Compatibility of mother wavelet functions with the electrocardiographic signal. 
                % In: . 2013 International Conference on Informatics, Electronics Vision (ICIEV) p. 1–4. 
                % http://dx.doi.org/10.1109/ICIEV.2013.6572707.    

            tic
            qmf=MakeONFilter(parameters.heart.reconstr_motherWavelet, parameters.heart.reconstr_Par);
            ECG_reconstr_A1 = recdecompsh(ECG,qmf);
            symmletReconstrTime = toc;

            if debugPlotON == 1           
                indices = [4 5 6]; colorPlot = [0.3 0.9 0.4]; 
                titleStr = ['Symmlet reconstruction, Par=', num2str(parameters.heart.reconstr_Par), ', t = ', num2str(symmletReconstrTime, '%3.1f'), ' sec'];
                pre_plotQuickECG(t, ECG_reconstr_A1, zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
                indices = [7 8 9]; colorPlot = [0.01 0.025 0.005]; 
                titleStr = ['DIFFERENCE: Input - Symmlet'];
                pre_plotQuickECG(t, ECG-ECG_reconstr_A1', zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
                drawnow
            end
        else
            disp(['           - Skipping Wavelet Reconstruction of ECG Signal'])
            
        end

    %% WAVELET DENOISE with SOFT THRESHOLDING (from WavDen)
    
        if parameters.heart.skip_waveletSoftThresholding ~= 1
            
            disp('           + Wavelet Soft Thresholding Denoising of ECG Signal')
            % Denoise using the Translation Invariant procedure with soft thresholding.
            tic
            ECG_sThr = recTI(ECG','S');
            waveleSoftThresholdTime = toc;

            if debugPlotON == 1           
                indices = [10 11 12]; colorPlot = [0.88 0.3 0.55]; 
                titleStr = ['Wavelet Soft Thresholding, ', 't = ', num2str(waveleSoftThresholdTime, '%3.1f'), ' sec'];
                pre_plotQuickECG(t, ECG_sThr, zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
                indices = [13 14 15]; colorPlot = [0.01 0.025 0.005];
                titleStr = ['DIFFERENCE: Input - Soft Thresholding'];
                pre_plotQuickECG(t, ECG-ECG_sThr', zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
                drawnow
            end              
            
        else
            disp('           - Skipping Wavelet Soft Thresholding Denoising of ECG Signal')
        end
        
    %% PLOT
    
        if debugPlotON == 1           
            % Auto-SAVE
            try
                if handles.figureOut.ON == 1                     
                    dateStr = plot_getDateString(); % get current date as string          
                    fileNameOut = sprintf('%s%s', 'debug_ECG-Denoising_', strrep(handles.inputFile, '.bdf', ''), '.png');
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
        
    try 
        ECG = ECG_reconstr_A1;
    catch err
        err
        error('add full-proof if-else things as now I just output the reconstructionas the others are too time-consuming')
    end
    
    
    function pre_plotQuickECG(t, ECG, zoomWindow, Fs, colorPlot, titleStr, parameters, handles, rows, cols, indices)
      
        % original timeseries
        subplot(rows,cols,indices(1));
        plot(t, ECG, 'Color', colorPlot); 
            xlim([min(t) max(t)]); xlabel('Time [s]'); ylabel('\muV'); title(titleStr);
            ylim([min(ECG) max(ECG)])
        
        % zoomed timeseries
        subplot(rows,cols,indices(2));
        plot(t, ECG, 'Color', colorPlot); 
            xlim(zoomWindow); xlabel('Time [s]'); ylabel('\muV'); title('Zoom');
            ylim([min(ECG) max(ECG)])
        
        % frequency domain
        subplot(rows,cols,indices(3));
        
            % FFT
            freqRange = 0 : 0.1 : parameters.filter.bandPass_hiFreq;
            nOverlap = parameters.powerAnalysis.nOverlap; % in %  
            
            [f, ECGStruct.amplit, ECGStruct.PSD, ECGStruct.amplit_matrix, ECGStruct.PSD_matrix] = ...
                    analyze_powerSpectrum(ECG, Fs, 'Tukey', 0.10, Fs*10, Fs*10, freqRange, nOverlap, 'ECG');
        
            % Plot
            semilogy(f, (ECGStruct.PSD.mean), 'Color', colorPlot); 
            xlim([0 80]); xlabel('Freq [Hz]'); ylabel('Amplitude'); title('FFT');
            