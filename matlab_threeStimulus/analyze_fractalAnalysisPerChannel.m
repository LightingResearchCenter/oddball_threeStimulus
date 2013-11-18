function fractalAnalysis = analyze_fractalAnalysisPerChannel(EEGchannel, ch, Fs, parameters, handles)

    debugMatFileName = 'temp_EEG_MFDFA.mat';
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

    parameters.fractalEEG.downsampleFactor = 256 / 4096;    
    parameters.fractalEEG.nufft_mode = 'fast'; % 'direct' or 'fast'
    
    downSampl = [256];
    % downSampl = [1]
    
    % You can test the effect of downsampling to your analysis, but
    % typically you only want to use one value
    for i = 1 : length(downSampl)
        
        parameters.fractalEEG.downsampleFactor = downSampl(i) / 4096;
        
        nanIndices = isnan(EEGchannel);
        numberOfNans = sum(nanIndices);

        debugON = 1;

        L = length(EEGchannel);
        x = (0 : 1 : L-1)';

        % remove NaNs
        x_nonNans = x(~nanIndices);
        EEG_nonNans = EEGchannel(~nanIndices);

        %% Now we have a problem that we have non-uniformly sampled EEG data
        % with the missing values              
        
            % alternative to downsampling and interpolation, you could try
            % to find artifact-free epochs
        
            % downsample    
            Fs_new = Fs * parameters.fractalEEG.downsampleFactor;
            newLength = L * parameters.fractalEEG.downsampleFactor;    
            xi = (0 : 1/parameters.fractalEEG.downsampleFactor : L-1)';

                % NUFFT
                % one way to go around the missing values is to use NUFFT reconstruction as for
                % example provided for Matlab by CMCL: http://www.cims.nyu.edu/cmcl/nufft/nufft.html
                %{
                nj = length(EEG_nonNans); % number of sources   (integer)
                ms = round(2.1*nj); % number of Fourier modes computed (-ms/2 to (ms-1)/2 )
                iflag = +1; % determines sign of FFT (see above)

                xj = x_nonNans; % location of sources (real *8) on interval [-pi,pi].
                cj = EEG_nonNans; % strengths of sources (complex *16)

                eps=1e-12;

                if strcmp(parameters.fractalEEG.nufft_mode, 'direct')
                    fk = dirft1d1(nj,xj,cj,iflag,ms); % "accurate" approximation            
                elseif strcmp(parameters.fractalEEG.nufft_mode, 'fast')    
                    fk = nufft1d1(nj,xj,cj,iflag,eps,ms); % faster approximation
                else
                    modeNow = parameters.fractalEEG.nufft_mode
                    error('Typo with the the NUFFT mode')
                end

                yi = ifft(fk);
                %}

                % DMD
                % http://userver.ftw.at/~vogel/NOS.html

                % Extended Fourier analysis of signals (EDFT)
                % http://arxiv.org/abs/1303.2033
                % http://www.mathworks.com/matlabcentral/fileexchange/11020-extended-dft

                % spline, quick'n'dirty
                EEG_downsampled = interp1(x_nonNans,EEG_nonNans,xi,'pchip');
                % [EEG_downsampled, xi] = pre_correctHeartRatePeriodForOutliers([], xi, EEG_downsampled, handles);
                EEG_downsampled = EEG_downsampled(~isnan(EEG_downsampled));
                xi = xi(~isnan(EEG_downsampled));


                % Nyquist interpolation
                % EEG_downsampled = interpft(EEGchannel, newLength);
        
                numberOfNans_new = sum(isnan(EEG_downsampled));

        %% PLOT INPUT
        
            if handles.figureOut.debugON == 1  

                scrsz = get(0,'ScreenSize'); % get screen size for plotting
        
                    fig = figure('Color', 'w');
                        set(fig, 'Position', [0.05*scrsz(3) 0.20*scrsz(4) 0.92*scrsz(3) 0.82*scrsz(4)])
                        sP(1) = subplot(2,2,[1 2]);
                            plot(xi/Fs, EEG_downsampled, 'r', x/Fs, EEGchannel, 'b'); 

                            title(['Number of NaNs from: "', num2str(numberOfNans), '" (', num2str(100*numberOfNans/L), '%) to "', num2str(numberOfNans_new), '"'])
                            legend(['Downsampled (@ ', num2str(Fs_new), ' Hz) and pchip interpolated'], 'Raw Input')
                            legend('boxoff')                    
                            xlabel('Time [s]')
                            xlim([min(x/Fs) max(x/Fs)])
                            drawnow

            end

        %% Do the MFDFA analyis
        
            disp(['              - Multifractal DFA (MFDFA)'])
            m = 1;        
            scmin = Fs_new;        
            scmax = newLength;
            ressc = 100;
            qmin = -5;
            qmax = 5;
            qres = abs(((qmin - qmax) / 0.2))+1;    
            % eps = 1; % 1 millisecond with HRV
            try
                tic                
                [s,q,Hq,h,Dh,logFq] = MFDFA(EEG_downsampled',m,scmin,scmax,ressc,qmin,qmax,qres);
                timing = toc;
            catch err            
                if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                    err
                    error('Have you added the MFDFA with subfolders to your Matlab path, as the function was not found now!')
                else
                    err
                    err.identifier
                    error('Some other error here, add cases while you get new errors!')
                end
            end

        %% Output
        
        % Variables of interest from Zorick and Mandelkern (2013)
        % http://dx.doi.org/10.1371/journal.pone.0068360

            % see short description of how these were calculated from the 2nd
            % column of page 2 of the Zorick and Mandelkern
            fractalAnalysis.scalar.MFDFA_mean_h = mean(h);
            fractalAnalysis.scalar.MFDFA_mean_Dh = mean(Dh);
            fractalAnalysis.scalar.MFDFA_width_h = max(h) - min(h);
            fractalAnalysis.scalar.MFDFA_height_Dh = max(Dh) - min(Dh);

        %% PLOT IF YOU WANT THE DEBUG
        if handles.figureOut.debugON == 1                      

            % as Figure 4 and Figure 5 of Zorick and Mandelkern (2013)
            % Zorick T, Mandelkern MA. 2013. 
            % Multifractal Detrended Fluctuation Analysis of Human EEG: Preliminary Investigation and Comparison with the Wavelet Transform Modulus Maxima Technique. 
            % PLoS ONE 8:e68360. 
            % http://dx.doi.org/10.1371/journal.pone.0068360

            sP(2) = subplot(2,2,3);
            p = plot(h, Dh , 'ko');
            xlabel('h (Hölder/Hurst Exponent)')
            ylabel('D(h) (FractalDimension/MultifractalSpectrum)')
            set(p, 'markerFaceColor', 'b')

            hold on
            l = line([fractalAnalysis.scalar.MFDFA_mean_h fractalAnalysis.scalar.MFDFA_mean_h], [0 1]);
            hold off

            xlims = get(gca, 'XLim');
            ylims = get(gca, 'YLim');
            yOffset = 0.1;

            tx(1) = text(0.95*xlims(2), ylims(2)*(1-yOffset*1), ['mean_h = ', num2str(fractalAnalysis.scalar.MFDFA_mean_h)]);
            tx(2) = text(0.95*xlims(2), ylims(2)*(1-yOffset*2), ['mean_D(h) = ', num2str(fractalAnalysis.scalar.MFDFA_mean_Dh)]);
            tx(3) = text(0.95*xlims(2), ylims(2)*(1-yOffset*3), ['width_h = ', num2str(fractalAnalysis.scalar.MFDFA_width_h)]);
            tx(4) = text(0.95*xlims(2), ylims(2)*(1-yOffset*4), ['height_D(h) = ', num2str(fractalAnalysis.scalar.MFDFA_height_Dh)]);

            tx(5) = text(0.95*xlims(2), ylims(2)*(1-yOffset*6), ['F_s = ', num2str(Fs_new), ' Hz (', num2str(newLength), ' samples)']);
            tx(6) = text(0.95*xlims(2), ylims(2)*(1-yOffset*7), ['t_{comp} = ', num2str(timing/60), ' min']);
            tx(7) = text(0.95*xlims(2), ylims(2)*(1-yOffset*8), ['scales = ', num2str(scmin), ':', num2str(scmax), ' [n = ', num2str(ressc), ']']);
            tx(8) = text(0.95*xlims(2), ylims(2)*(1-yOffset*9), ['m = ', num2str(m)]);
            set(tx, 'HorizontalAlignment', 'right')

            set(sP, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)
            set(tx, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)

             % Auto-SAVE
            handles.figureOut.resolution = '-r150';
            try
                if handles.figureOut.ON == 1      
                    drawnow
                    dateStr = plot_getDateString(); % get current date as string          
                    fileNameOut = sprintf('%s%s', 'MFDFA_v2_EEG_ch', num2str(ch), '_m', num2str(m), '_', num2str(Fs_new), 'Hz_', num2str(newLength), 'samples_', strrep(handles.inputFile, '.bdf', ''), '.png');
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
    end