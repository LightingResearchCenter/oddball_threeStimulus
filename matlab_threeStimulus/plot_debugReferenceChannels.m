function plot_debugReferenceChannels(ref, refRightEar, refLeftEar, parameters, handles)
        
    style = handles.style;        
    x = (1 : 1 : length(ref))';

    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color','w');      
        set(fig, 'Position', [0.10*scrsz(3) 0.05*scrsz(4) 0.82*scrsz(3) 0.9*scrsz(4)])

    rows = 3;
    cols = 2;

    %% Input Time-domain
    i = 1;
    s(i) = subplot(rows,cols,i);
        p(i,:) = plot(x, ref, 'k', x, refRightEar, 'r', x, refLeftEar, 'b');
        xlabel('Samples'); ylabel('\muV');    
        xlim([min(x) max(x)])
        leg = legend('Ref', 'Ref rightEar', 'Ref leftEar', 3, 'Location', 'Best');
            legend('boxoff')
            tit(i) = title('Unfiltered Reference Input');
            drawnow

    %% Compute power spectrum

        % set defaults, and trim the input
        Fs = parameters.EEG.srate;
        [~, ~, ~, ~, ~, ~, ref] = debug_setPowerAnalysisDefaults(ref, Fs, parameters);
        [~, ~, ~, ~, ~, ~, refRightEar] = debug_setPowerAnalysisDefaults(refRightEar, Fs, parameters);
        [windowType, r, segmentLength, nfft, freqRange, nOverlap, refLeftEar] = debug_setPowerAnalysisDefaults(refLeftEar, Fs, parameters);

        % compute for each three vectors
        [~, refStruct.amplit, refStruct.PSD, refStruct.amplit_matrix, refStruct.PSD_matrix] = ...
                analyze_powerSpectrum(ref, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');
        [~, refRightEarStruct.amplit, refRightEarStruct.PSD, refRightEarStruct.amplit_matrix, refRightEarStruct.PSD_matrix] = ...
                analyze_powerSpectrum(refRightEar, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');
        [f, refLeftEarStruct.amplit, refLeftEarStruct.PSD, refLeftEarStruct.amplit_matrix, refLeftEarStruct.PSD_matrix] = ...
                analyze_powerSpectrum(refLeftEar, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');


    %% Plot the spectrum analysis
    i = 2;
    s(i) = subplot(rows,cols,i);

        hold on
        %{
        p(i,1) = plot(f, refStruct.amplit.mean, 'k');
        p(i,2) = plot(f, refRightEarStruct.amplit.mean, 'r');
        p(i,3) = plot(f, refLeftEarStruct.amplit.mean, 'b');           
        %}
        p(i,1) = plot(f, 10*log10(refStruct.PSD.mean), 'k');
        p(i,2) = plot(f, 10*log10(refRightEarStruct.PSD.mean), 'r');
        p(i,3) = plot(f, 10*log10(refLeftEarStruct.PSD.mean), 'b');                       
        hold off

        xlabel('Hz'); ylabel('dB');                        
        xlim([0 70])
        % ylim([0 20])
        tit(i) = title('PSD of Unfiltered Input');
        drawnow

    %% Apply a notch filter (from BioSig)

        % Mode    
        %		'PCA 60'	60 Hz PCA/SVD filter
        %		'NOTCH 60' 	60 Hz FIR Notch filter, order=3
        %       'FIT 60'        fit and remove 60 Hz sine/cosine wave
        HDR.SampleRate = Fs;
        ref_pca60 = pre_remove5060hz_modified(ref, HDR, 'PCA 60');
        ref_pca60_x2 = pre_remove5060hz_modified(ref_pca60, HDR, 'PCA 60'); % filter twice
        ref_notch60 = pre_remove5060hz_modified(ref, HDR, 'NOTCH 60');
        ref_fit60 = pre_remove5060hz_modified(ref, HDR, 'FIT 60');


    %% Power analysis of the filtered signal

        % compute for each three vectors
        [~, PCA60.amplit, PCA60.PSD, PCA60.amplit_matrix, PCA60.PSD_matrix] = ...
                analyze_powerSpectrum(ref_pca60, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');
        [~, PCA60_x2.amplit, PCA60_x2.PSD, PCA60_x2.amplit_matrix, PCA60_x2.PSD_matrix] = ...
                analyze_powerSpectrum(ref_pca60_x2, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');
        [~, notch60.amplit, notch60.PSD, notch60.amplit_matrix, notch60.PSD_matrix] = ...
                analyze_powerSpectrum(ref_notch60, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');
        [f, fit60.amplit, fit60.PSD, fit60.amplit_matrix, fit60.PSD_matrix] = ...
                analyze_powerSpectrum(ref_fit60, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');


    %% Plot the filtered PSD
    i = 3;
    s(i) = subplot(rows,cols,i);

        hold on
        %{
        p(i,1) = plot(f, refStruct.amplit.mean, 'k');
        p(i,2) = plot(f, refRightEarStruct.amplit.mean, 'r');
        p(i,3) = plot(f, refLeftEarStruct.amplit.mean, 'b');           
        %}
        p(i,1) = plot(f, 10*log10(PCA60.PSD.mean), 'g');
        p(i,2) = plot(f, 10*log10(PCA60_x2.PSD.mean), 'k--');
        p(i,3) = plot(f, 10*log10(notch60.PSD.mean), 'm');
        p(i,4) = plot(f, 10*log10(fit60.PSD.mean), 'Color', [0 0.749 0.749]);                       
        hold off

        xlabel('Hz'); ylabel('dB');                        
        xlim([0 70])
        % ylim([0 20])
        tit(i) = title('PSD of Notch-filtered Input (Ref.)');
            leg = legend('PCA 60', 'PCA 60 x2', 'NOTCH 60', 'FIT 60', 4, 'Location', 'Best');
                    legend('boxoff')
        drawnow



    %% Plot the filtered time-domain
    x2 = (1 : 1 : length(ref))';
    whos

    i = 4;
    s(i) = subplot(rows,cols,i);
        p(i,:) = plot(x2, ref_pca60, 'g', x2, ref_pca60_x2, 'k--', x2, ref_notch60, 'y', x2, ref_fit60, 'k');
            set(p(i,4), 'Color', [0 0.749 0.749])
            xlabel('Samples'); ylabel('\muV');    
            xlim([min(x2) max(x2)])                
                titstr = sprintf('%s\n%s', '60Hz removed Reference', 'with three methods from BioSig');
                tit(i) = title(titstr);
                drawnow

    %% Filter with band-pass (GENERAL)
    loFreq = parameters.filter.bandPass_loFreq;
    hiFreq = parameters.filter.bandPass_hiFreq;
    filterOrder = parameters.filterOrder;
    ref_bpF = pre_bandbassFilter(ref, parameters.EEG.srate, [hiFreq, loFreq], filterOrder);
    ref_pca60_bp = pre_bandbassFilter(ref_pca60, parameters.EEG.srate, [hiFreq, loFreq], filterOrder);

    %% Power analysis of the Bandpass-FILTERED

        [~, PCA60_bp.amplit, PCA60_bp.PSD, PCA60_bp.amplit_matrix, PCA60_bp.PSD_matrix] = ...
                analyze_powerSpectrum(ref_pca60_bp, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');
        [f, ref_bp.amplit, ref_bp.PSD, ref_bp.amplit_matrix, ref_bp.PSD_matrix] = ...
                analyze_powerSpectrum(ref_bpF, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugReference');

    %% PLOT POWER SPECTRUM
    i = 5;
    s(i) = subplot(rows,cols,i);

        hold on
        %{
        p(i,1) = plot(f, refStruct.amplit.mean, 'k');
        p(i,2) = plot(f, refRightEarStruct.amplit.mean, 'r');
        p(i,3) = plot(f, refLeftEarStruct.amplit.mean, 'b');           
        %}
        p3(i,1) = plot(f, 10*log10(PCA60_bp.PSD.mean), 'k');
        p3(i,2) = plot(f, 10*log10(ref_bp.PSD.mean), 'r');                        
        p3(i,3) = plot(f, 10*log10(refStruct.PSD.mean), 'b');  
        hold off

        xlabel('Hz'); ylabel('dB');                        
        xlim([0 70])
        % ylim([0 20])
        tit(i) = title(['PSD of Bandpass-filtered Signals [', num2str(loFreq), '-', num2str(hiFreq), ' Hz], order: ', num2str(filterOrder)]);
            leg = legend('PCA 60', 'Ref. BandPass', 'Ref', 3, 'Location', 'Best');
                    legend('boxoff')
        drawnow



    %% Plot the filtered time-domain        
    whos

    i = 6;
    s(i) = subplot(rows,cols,i); 
        hold on
        p3(i,1) = plot(x2, ref_bpF, 'r');
        p3(i,2) = plot(x2, ref, 'b');
        p3(i,3) = plot(x2, ref_pca60_bp, 'k');
        hold off
            xlabel('Samples'); ylabel('\muV');    
            xlim([min(x2) max(x2)])                
            %xlim([254000 255000])
            %ylim([-55 55])
                titstr = sprintf('%s\n%s', 'Time series of Bandpass-filtered Signals', ['[', num2str(loFreq), '-', num2str(hiFreq), ' Hz], order: ', num2str(filterOrder)]);
                tit(i) = title(titstr);
                drawnow

    %% General Styling
    set(s, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
    set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
    set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase+1, 'FontWeight', 'bold')
    drawnow

    %% Auto-SAVE
    try
        if handles.figureOut.ON == 1                      
            dateStr = plot_getDateString(); % get current date as string          
            fileNameOut = sprintf('%s%s', 'plot_refRemoval_', strrep(handles.inputFile, '.bdf', ''), '.png');
            export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            %cd(path.code)
        end
    catch err
        err
        str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
        error(str)
    end

    pause