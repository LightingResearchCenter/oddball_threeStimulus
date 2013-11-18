function basic_cleaningOfRecording()

    % Just band-pass filters the data and corrects for EOG/ECG artifacts
    
        % Useful if you want to analyze the cleaned signal with some other
        % software

    tic

    % Petteri Teikari, petteri.teikari@gmail.com, 2013
    % Lighting Research Center, Rensselaer Polytechnic Institute, Troy, NY
    close all
    clear all
    
    %% General Settings
    
        % i.e. like where the folders are, fonts to be used, colors, etc.
        handles = init_DefaultSettings(); % use a subfunction        
        
        plotPhases = 1;
        scrsz = handles.style.scrsz;
        fig = figure('Color', 'w');
            set(fig, 'Position', [0.1*scrsz(3) 0.05*scrsz(4) 0.65*scrsz(3) 0.90*scrsz(4)])
            rows = 5;
            cols = 3;
            zoomX = 0.4; % fraction of all x samples
            xDuration = 0.5; % seconds
        

    %% Parameters for ANALYSIS
    
        % i.e. like artifact rejection thresholds, filter cutoffs,
        % numerical parameters for EEG analysis
        handles.parameters = init_DefaultParameters(handles); % use a subfunction            
        
    %% INPUT
    
        filesIn = {'da_01_10.bdf'};
        fileNameIn = filesIn;        
        handles.inputFile = fileNameIn;
        inputFiles = fullfile(handles.path.dataFolder, fileNameIn);
        if ~iscell(inputFiles) % for only one file
            numberOfFiles = 1;
        else
            numberOfFiles = length(inputFiles);
        end
        i = 1;
        
    %% IMPORT THE DATA
     
        % You could modify what is inside here if you have EEG recorded
        % with some other system than with BioSemi ActiveTwo
        handles.inputFile = filesIn{i}; 
        inputFiles = fullfile(handles.path.dataFolder, filesIn{i});
        [dataMatrix, triggers, info, handles.parameters.EEG.srate] = IMPORT_eegData(i, fileNameIn, inputFiles, handles);        
        
        t = linspace(1, length(dataMatrix), length(dataMatrix)) / handles.parameters.EEG.srate;

    %% CONDITION THE DATA    
    
        [dataSamples, colsIn] = size(dataMatrix);
        HDR.SampleRate = handles.parameters.EEG.srate;
        
            if plotPhases == 1            
                i = 1;
                [p, sp] = plotRow([], [], rows, cols, i, t, dataMatrix(:, 1+handles.parameters.BioSemi.chOffset:handles.parameters.EEG.nrOfChannels+handles.parameters.BioSemi.chOffset), handles, zoomX, xDuration, 'Raw Input');
                drawnow                    
            end
    
        for j = 1 : colsIn     
            
            disp([num2str(j), ':th ', 'channel'])
            
            disp(['... notch filter']);
            dataMatrix(:,j) = pre_remove5060hz_modified(dataMatrix(:,j), HDR, 'PCA 60');                
            
            disp(['... bandpass filter (', num2str(handles.parameters.filter.bandPass_loFreq), '-', num2str(handles.parameters.filter.bandPass_hiFreq), ' Hz), order = ', num2str(handles.parameters.filterOrder)]);
            dataMatrix(:,j) = pre_bandbassFilter(dataMatrix(:,j), handles.parameters.EEG.srate, [handles.parameters.filter.bandPass_hiFreq handles.parameters.filter.bandPass_loFreq], handles.parameters.filterOrder, 0, handles);
            
        end
        
            if plotPhases == 1            
                i = 4;
                [p, sp] = plotRow(p, sp, rows, cols, i, t, dataMatrix(:, 1+handles.parameters.BioSemi.chOffset:handles.parameters.EEG.nrOfChannels+handles.parameters.BioSemi.chOffset), handles, zoomX, xDuration, 'Filtered');
                drawnow                    
            end
                
            
        disp(' remove reference (the ears)')
        REF = dataMatrix(:, (1:handles.parameters.BioSemi.chOffset));
        EEG = dataMatrix(:, (1+handles.parameters.BioSemi.chOffset:handles.parameters.EEG.nrOfChannels+handles.parameters.BioSemi.chOffset));
        EOG = dataMatrix(:,7);
        ECG = dataMatrix(:,8);            
        reference = (REF(:,1) + REF(:,2) ) / 2;
        for j = 1 : handles.parameters.EEG.nrOfChannels
            EEG(:,j) = EEG(:,j) -reference;
        end
        EOG = EOG -reference;
        ECG = ECG -reference;            
        
            if plotPhases == 1            
                i = 7;
                [p, sp] = plotRow(p, sp, rows, cols, i, t, EEG, handles, zoomX, xDuration, 'referenceRemoved');
                drawnow                    
            end

        % EOG CORRECTION                   
        
            disp('  correct for EOG artifacts')

            % combine EOG and EEG
            S1 = [EEG EOG];
            EL = [1 2 3 4];            
            OL = 5;
            try
                R = regress_eog(S1, EL, OL);
            catch err
                err
                if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                    error('EOG Rejection not done! Have you installed the BioSig toolbox needed for EOG artifact rejection? http://sourceforge.net/projects/biosig/')
                else
                    err
                end                
            end

            % apply the correction, with offset correction
            S2 = [ones(size(S1,1),1),S1] * R.r1;
            % S2b = S1 * R.r0;    % without offset correction
            EEG = S2(:,EL);                     
            
            if plotPhases == 1            
                i = 10;
                [p, sp] = plotRow(p, sp, rows, cols, i, t, EEG, handles, zoomX, xDuration, 'EOG Corrected');
                drawnow                    
            end


        % ECG CORRECTION

            disp('   correct for ECG artifacts')
            S1 = [EEG ECG];
            EL = [1 2 3 4];
            OL = 5;
            R = regress_eog(S1, EL, OL);
            S2 = [ones(size(S1,1),1),S1] * R.r1;
            EEG = S2(:,EL);
            
            if plotPhases == 1            
                i = 13;
                [p, sp] = plotRow(p, sp, rows, cols, i, t, EEG, handles, zoomX, xDuration, 'ECG Corrected');
                drawnow                    
            end

      
    %% OUTPUT
    
        % assign to more intuitive variables
        fileNameOut = sprintf('%s', strrep(filesIn{1}, '.bdf', '_cleaned.mat'));
        buttonTrigger = triggers.button;
        whos
        save(fullfile(handles.path.dataFolder, fileNameOut), 'EEG', 'EOG', 'ECG', 'reference', 'buttonTrigger')
        
        
    %% DOWNSAMPLE if needed
    
    
    %% Auto-SAVE

        try
            if handles.figureOut.ON == 1      
                drawnow
                dateStr = plot_getDateString(); % get current date as string        
                
                loString = num2str(handles.parameters.filter.bandPass_loFreq);
                
                if handles.parameters.filter.bandPass_hiFreq < 10
                    hiString = ['0', num2str(handles.parameters.filter.bandPass_hiFreq)];
                else
                    hiString = num2str(handles.parameters.filter.bandPass_hiFreq);
                end
                
                fileNameOut = sprintf('%s%s%s%s %s%s%s%s %s%s%s%s %s', 'cleaningPlot_', ...                    
                    '_f0Lo-', loString, 'Hz', ...
                    '_f0Hi-', hiString, 'Hz', ...                                        
                    '_', strrep(handles.inputFile, '.bdf', ''), '.png');
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
            end
        catch err
            err
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        end
    
    
    function [p, sp] = plotRow(p, sp, rows, cols, i, t, EEG, handles, zoomX, xDuration, titStr)
                
        sp(i) = subplot(rows,cols,i);
            p(i,:) = plot(t, EEG);
            tit = title(titStr);
             ylab(i) = ylabel('Amplitude [\muV]');
             xlab(i) = xlabel('Time [s]');
             xlim([1 t(end)])

        i = i+1;
        sp(i) = subplot(rows,cols,i);
            p(i,:) = plot(t, EEG);            
            XStart = t(round(zoomX*length(t)));
            xlim([XStart XStart+xDuration])                
            ylab(i) = ylabel('Amplitude [\muV]');
            xlab(i) = xlabel('Time [s]');
            
        % POWER SPECTRUM 
            
            Fs = handles.parameters.EEG.srate;
            segmentLength = handles.parameters.powerAnalysis.segmentLength * Fs;
            % segmentLength = segmentLength / 10;
            nfft = segmentLength;
            nOverlap = handles.parameters.powerAnalysis.nOverlap; % in %
            nOverlap = 50;
            r = handles.parameters.powerAnalysis.tukeyWindowR;       
            windowType = 'Tukey';
            freqRange = 0 : 0.1 : handles.parameters.filter.bandPass_hiFreq;

            % trim the channel data to be an integer of sampling frequency
            P = nextpow2(length(EEG));
            numberOfSegmentLengths = floor(length(EEG) / segmentLength);
            
            [noOfSamples, noOfChannels] = size(EEG);       
            for ch = 1 : noOfChannels
                EEGchannel = EEG(1:(numberOfSegmentLengths*segmentLength), ch);
                [f, amplit(:,ch), PSD, amplit_matrix, PSD_matrix] = analyze_powerSpectrum(EEGchannel, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'eegTimeSeries');                
            end      
            
            % put to a matrix the structure
            ps = zeros(length(amplit(:,1).mean), noOfChannels);
            for ch = 1 : noOfChannels
               ps(:, ch) = amplit(:,ch).mean;
            end
            
            
        % Plot Power spectrum
        i = i+1;
        sp(i) = subplot(rows,cols,i);            
        p(i,:) = plot(f, 10*log10(ps));    
        
            xlim([freqRange(1) 20])       
            ylab(i) = ylabel('Amplitude [dB]');
            xlab(i) = xlabel('Frequency [Hz]');
            tit2 = title('FFT');
            
            leg = legend('Fz', 'Cz', 'Pz', 'Oz', 'Location', 'NorthEastOutside');
                legend('boxoff')
                
            if strcmp(titStr, 'Filtered')
                yLims = get(gca, 'YLim');
                xLims = get(gca, 'XLim');
                txt = ['bandpass: ', num2str(handles.parameters.filter.bandPass_loFreq), '-', num2str(handles.parameters.filter.bandPass_hiFreq), ' Hz'];
                xRatio = 0.97; yRatio = 0.90;
                tx = text(xLims(2)*xRatio, yLims(2)*yRatio, txt, 'HorizontalAlignment', 'right');
                set(tx, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)
            end
                    
                
        % style
        set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)
        set([tit tit2], 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase+1, 'FontWeight', 'bold')
        %set(ylab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
            
            
        