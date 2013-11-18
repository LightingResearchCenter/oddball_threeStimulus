function plot_debugBandPassPlot(data, dataOut, dataOut2, ferror, b, a, v, frec3dB_final, xdB_at_fx, orderx, Fs, lowPass, highPass, N, typef, parameters, handles)

    rows = 2;
    cols = 2;
    style = handles.style;

    close all

    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color','w');      
        set(fig, 'Position', [0.10*scrsz(3) 0.05*scrsz(4) 0.52*scrsz(3) 0.6*scrsz(4)])

    
    % do spectral analysis
    [~, ~, ~, ~, ~, ~, data] = debug_setPowerAnalysisDefaults(data, Fs, parameters);
    [windowType, r, segmentLength, nfft, freqRange, nOverlap, dataOut] = debug_setPowerAnalysisDefaults(dataOut, Fs, parameters);

    % compute for both input and output signals
    [~, dataInPS.amplit, dataInPS.PSD, dataInPS.amplit_matrix, dataInPS.PSD_matrix] = ...
            analyze_powerSpectrum(data, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugBandpass');
    [f, dataOutPS.amplit, dataOutPS.PSD, dataOutPS.amplit_matrix, dataOutPS.PSD_matrix] = ...
            analyze_powerSpectrum(dataOut, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, 'debugBandpass');

    x = (linspace(1,length(data),length(data)))';
    t = x / Fs;

        ind = 1;
        sp(ind) = subplot(rows,cols,ind);
            hold on
            plot(t, data, 'r')            
            plot(t, dataOut, 'k')
            hold off
            xlabel('Time [s]')
            ylabel('\muV')
            tit(ind) = title('Time Domain');
            leg(ind) = legend('Input', 'Output', 2, 'Location', 'Best');
                legend('boxoff')

        ind = 2;
        sp(ind) = subplot(rows,cols,ind);
            hold on
            plot(f, 10*log10(dataInPS.PSD.mean), 'r')            
            plot(f, 10*log10(dataOutPS.PSD.mean), 'k')
            xlabel('Hz')
            ylabel('dB')
            ylim([-50 50])
            hold off
            grid on
            titStr = sprintf('%s\n%s', 'Power Spectral Density', ['ERPLab, Zero-Phase IIR Butterworth, order ', num2str(N), ' (', num2str(lowPass), '-', num2str(highPass), ' Hz)']);
            tit(ind) = title(titStr);
            leg(ind) = legend('Input', 'Output', 2, 'Location', 'Best');
                legend('boxoff')


        freqMax = 70;  
        % The frequency/phase response of the filter

        n = logspace(-10,1);
        [h_lo, w_lo] = freqz(b(1,:),a(1,:), n);        
        [h_hi, w_hi] = freqz(b(2,:),a(2,:), n);

        mag_lo = 20*log10(abs(h_lo));
        phase_lo = angle(h_lo);
        mag_hi = 20*log10(abs(h_hi));
        phase_hi = angle(h_hi);
        
        % combine hi and lo
        mag = mag_lo - mag_hi;
        phase = phase_lo - phase_hi;
        
        % convert angular frequency, first to normalized frequency and then
        % weigh with the Nyquist frequency (sampling rate / 2)
        f = w_lo/(2*pi);                
        f = f * (Fs/2);
        

        % MAGNITUDE
        ind = 3;
        sp(ind) = subplot(rows,cols,ind);
            plot(f,mag_lo,'g',f,mag_hi,'b')
            grid on
            xlabel('Frequency [Hz]')
            ylabel('dB')
            titStr = sprintf('%s\%s', 'Magnitude response', [' ']);
            tit(ind) = title(titStr);
            leg(ind) = legend('LowPass', 'HighPass', 2, 'Location', 'Best');
                legend('boxoff')

        % PHASE
        ind = 4;
        sp(ind) = subplot(rows,cols,ind);
            plot(f,phase_lo,'g',f,mag_lo,'b')
            grid on
            xlabel('Frequency [Hz]')
            ylabel('degree')
            titStr = sprintf('%s\%s', 'Phase response', [' ']);
            tit(ind) = title(titStr);
            leg(ind) = legend('LowPass', 'HighPass', 2, 'Location', 'Best');
                legend('boxoff')

        set(sp(1), 'XLim', [min(t) max(t)])
        set(sp(2:end), 'XLim', [0 freqMax])

        set(sp, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
        set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
        set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase+1, 'FontWeight', 'bold')

         % Auto-SAVE
        try
            if handles.figureOut.ON == 1      

                % correct number strings
                if N < 10
                    orderStr = ['00', num2str(N)];
                elseif N < 100
                    orderStr = ['0', num2str(N)];
                else
                    orderStr = num2str(N);                    
                end

                if lowPass < 10
                    lowStr = ['0', num2str(lowPass)];
                else
                    lowStr = num2str(lowPass);
                end

                drawnow
                dateStr = plot_getDateString(); % get current date as string          
                fileNameOut = sprintf('%s%s', 'plot_BandPassBehavior_order', orderStr, '_freqRange_', lowStr, '-', num2str(highPass), 'Hz_', strrep(handles.inputFile, '.bdf', ''), '.png');
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
            end
        catch err
            err
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        end