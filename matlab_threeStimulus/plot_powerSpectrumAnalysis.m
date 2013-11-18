 function plot_powerSpectrumAnalysis(f, Amplit, PSD, Amplit_matrix, PSD_matrix, channels, style, parameters, handles)
        
    debugMatFileName = 'tempPowerSpecPlot.mat';
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

        segmentLength = parameters.powerAnalysis.segmentLength * parameters.EEG.srate;
        x = 1 : 1 : segmentLength/2 + 1;
        freqRange = f;
        chNameOffset = 2; % no ref channels


    %% 2D Mean Plots

        scrsz = get(0,'ScreenSize'); % get screen size for plotting
        fig = figure('Color','w');      
            set(fig, 'Position', [0.05*scrsz(3) 0.15*scrsz(4) 0.62*scrsz(3) 0.9*scrsz(4)])

        % subplot layout
        rows = 1 * length(Amplit);
        cols = 2;                


        for ch = 1 : length(Amplit)            

            % PSD                        
            ind = ((ch-1) * cols) + 1;
            sp(ind) = subplot(rows, cols, ind);
                hold on                    
                pSD(ind,:) = plot(f, 10*log10(PSD{ch}.mean-PSD{ch}.SD), 'r', f, 10*log10(PSD{ch}.mean+PSD{ch}.SD), 'r');
                pMean(ind) = plot(f, 10*log10(PSD{ch}.mean), 'b');
                hold off
                    xlim([min(freqRange) parameters.filter.bandPass_hiFreq])

                    titStr = sprintf('%s\n%s', ['Power Spectral Density - ', parameters.BioSemi.chName{ch+chNameOffset}],...
                        ['segment length = ', num2str(parameters.powerAnalysis.segmentLength), ' sec, Overlap = ', num2str(parameters.powerAnalysis.nOverlap), '%']);
                    tit(ind) = title(titStr);
                    xlabel('Hz'); ylabel('dB');
                    grid on;

            % Amplit
            ind = ((ch-1) * cols) + 2;
            sp(ind) = subplot(rows, cols, ind);
                hold on                    
                pSD(ind,:) = plot(f, (Amplit{ch}.mean-Amplit{ch}.SD), 'r', f, (Amplit{ch}.mean+Amplit{ch}.SD), 'r');
                pMean(ind) = plot(f, Amplit{ch}.mean, 'b');
                hold off
                    xlim([min(freqRange) parameters.filter.bandPass_hiFreq])

                    titStr = sprintf('%s\n%s', ['Amplitude Spectrum - ', parameters.BioSemi.chName{ch+chNameOffset}],...
                        ['segment length = ', num2str(parameters.powerAnalysis.segmentLength), ' sec, Overlap = ', num2str(parameters.powerAnalysis.nOverlap), '%']);
                    tit(ind) = title(titStr);
                    xlabel('Hz'); ylabel('dB');
                    grid on;

                    leg = legend([pMean(1), pSD(1,1)], 'Mean', 'SD');
                        legend('boxoff')

        end


        set(sp, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
        set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
        set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase+1, 'FontWeight', 'bold')

        set(pMean, 'LineWidth', 2, 'Color', [0 0 0])              
        set(pSD, 'LineWidth', 1, 'Color', [0.9888 0.486 0.612])

        % Auto-SAVE
        try
            if handles.figureOut.ON == 1      
                drawnow
                dateStr = plot_getDateString(); % get current date as string          
                fileNameOut = sprintf('%s%s', 'plot_powerSpec_segmL-', num2str(parameters.powerAnalysis.segmentLength), 'sec_', strrep(handles.inputFile, '.bdf', ''), '.png');
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
            end
        catch err
            err
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        end

    %% TIME-FREQUENCY CONTOUR PLOT
    %{
            fig = figure('Color','w');      
                set(fig, 'Position', [0.095*scrsz(3) 0.60*scrsz(4) 0.62*scrsz(3) 0.9*scrsz(4)])

                % subplot layout
                rows = 2; % 1 * length(Amplit);
                cols = 2;                      
                [dataSamples, nrOfSegments] = size(Amplit_matrix{1});

                segmentVector = 1 : 1 : nrOfSegments;

                % freq indices
                i1 = find((f == min(freqRange)) == 1)
                if isempty(i1)
                   i1 = 1; 
                end
                i2 = find((f == parameters.filter.bandPass_hiFreq) == 1)

                n = 128;
                colorBarLimits_Amplit = [0 4.5]; 
                colorBarLimits_PSD = [-10 10]; 

                for ch = 1 : 2 % length(Amplit)

                    % PSD                        
                    ind = ((ch-1) * cols) + 1;
                    sp(ind) = subplot(rows, cols, ind);            

                        axis on                    

                        [c,h(ind)] = contourf(f(i1:i2), segmentVector, 10*log10(PSD_matrix{ch}(i1:i2,:))', n);                  

                        xlim([min(freqRange) parameters.filter.bandPass_hiFreq])
                        ylim([min(segmentVector) max(segmentVector)])

                        xlab = xlabel('Hz');
                        zlab = ylabel('Segments');
                        ylab = zlabel('dB');                    

                         titStr = sprintf('%s\n%s', ['Power Spectral Density - ', parameters.BioSemi.chName{ch+chNameOffset}],...
                            ['segment length = ', num2str(parameters.powerAnalysis.segmentLength), ' sec, Overlap = ', num2str(parameters.powerAnalysis.nOverlap), '%']);
                        tit(ind) = title(titStr);

                        colorbar
                        caxis(colorBarLimits_PSD)



                    % Amplit
                    ind = ((ch-1) * cols) + 2;
                    sp(ind) = subplot(rows, cols, ind);

                        axis on                    

                        [c,h(ind)] = contourf(f(i1:i2), segmentVector,  Amplit_matrix{ch}(i1:i2,:)', n);                  

                        xlim([min(freqRange) parameters.filter.bandPass_hiFreq])
                        ylim([min(segmentVector) max(segmentVector)])

                        xlab = xlabel('Hz');
                        ylab = ylabel('Segments');
                        zlab = zlabel('dB');        

                        titStr = sprintf('%s\n%s', ['Amplitude Spectrum - ', parameters.BioSemi.chName{ch+chNameOffset}],...
                            ['segment length = ', num2str(parameters.powerAnalysis.segmentLength), ' sec, Overlap = ', num2str(parameters.powerAnalysis.nOverlap), '%']);
                        tit(ind) = title(titStr);

                        colorbar
                        caxis(colorBarLimits_Amplit)


                end

                set(sp, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
                set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase+1, 'FontWeight', 'bold')                    
                set(h,'Edgecolor','none') 

                % Auto-SAVE
                try
                    if handles.figureOut.ON == 1      
                        drawnow
                        dateStr = plot_getDateString(); % get current date as string          
                        fileNameOut = sprintf('%s%s', 'plot_powerSpec_contour_segmL-', num2str(parameters.powerAnalysis.segmentLength), 'sec_', strrep(handles.inputFile, '.bdf', ''), '.png');
                        export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                        %cd(path.code)
                    end
                catch err
                    err
                    str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                                  'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
                    error(str)
                end

    %}

    %% NOISE ESTIMATE PLOT


        if strcmp(handles.inputFile, 'petteri_testSession6cycles_wECG.bdf')
            fig = figure('Color','w');      
                set(fig, 'Position', [0.075*scrsz(3) 0.60*scrsz(4) 0.62*scrsz(3) 0.225*scrsz(4)])

                % For the given debug session, the electrodes for Pz and Oz
                % were not attached to the scalp, and we can estimate noise
                % spectra from those
                rows = 1;
                cols = 2;

                % PSD    
                ind = 1;
                sp(ind) = subplot(rows, cols, ind);

                    noiseEstimate = (10*log10(PSD{3}.mean) + 10*log10(PSD{4}.mean)) / 2;
                    p2(ind,:) = plot(f, (10*log10(PSD{1}.mean))-noiseEstimate, 'b', f, (10*log10(PSD{2}.mean))-noiseEstimate, 'g');
                    xlim([min(freqRange) parameters.filter.bandPass_hiFreq])

                    tit(ind) = title(['Power Spectral Density']);
                    xlabel('Hz'); ylabel('dB');
                    grid on;

                    leg(ind) = legend('Fz', 'Cz');
                        legend('boxoff')

                % Amplit
                ind = 2;
                sp(ind) = subplot(rows, cols, ind);

                    noiseEstimate = (Amplit{3}.mean + Amplit{4}.mean) / 2;
                    p2(ind,:) = plot(f, Amplit{1}.mean-noiseEstimate, 'b', f, Amplit{2}.mean-noiseEstimate, 'g');
                    xlim([min(freqRange) parameters.filter.bandPass_hiFreq])

                    tit(ind) = title(['Amplitude Spectrum']);
                    xlabel('Hz'); ylabel('dB');
                    grid on;

                    leg(ind) = legend('Fz', 'Cz');
                        legend('boxoff')

                    set(sp, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
                    set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase, 'Location', 'Best')
                    set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase+1, 'FontWeight', 'bold')

                % Auto-SAVE
                    try
                        if handles.figureOut.ON == 1      
                            drawnow
                            dateStr = plot_getDateString(); % get current date as string          
                            fileNameOut = sprintf('%s%s', 'plot_powerSpec_noiseEstimatesegmL-', num2str(parameters.powerAnalysis.segmentLength), 'sec_', strrep(handles.inputFile, '.bdf', ''), '.png');
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


