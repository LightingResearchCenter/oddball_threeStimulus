function y = pre_bandbassFilter(data, Fs, cutOffs, filterOrder, filterOrderSteep, handles)

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempPreBandbassFilter.mat';
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
   
    
    %% FieldTrip
        %{
        dataMatrix_filt(:,i) = ft_preproc_bandpassfilter(dataMatrix(:,i), srate, ...
                            [handles.parameters.filter.bandPass_hiFreq, handles.parameters.filter.bandPass_loFreq], N, 'fir');
        %}

    %% EEGLAB
        %{
        % Requires Matlab Signal Processing Toolbox
        dataMatrix_filt(i,:) = eegfilt(dataMatrix(:,i)', srate, ...
                            [handles.parameters.filter.bandPass_loFreq, handles.parameters.filter.bandPass_hiFreq]);

        %}

    %% Old Code
    
        %{
        % Redefine the cutoff frequencies with the Nyquist criterion
        nyq = Fs/2;
        low2 = cutOffs(2) / nyq;
        high2 = cutOffs(1) / nyq;
        passband = [low2 high2];
        
        % Define filter function
        [Bb Ab] = butter(filterOrder, passband);

        % Actually filter
        filt = filter(Bb,Ab,data);
        %}
    
    %% ZERO-PHASE FILTER to preserve time domain characteristics
    
        % for short intro see:
        % http://erpinfo.org/erplab/erplab-documentation/documentation-archive-for-previous-versions/v1.x-documentation/erplab-manual/Filtering.htm
        
        % For more discussion see e.g.
        
            % Widmann A, Schröger E. 2012. 
            % Filter effects and filter artifacts in the analysis of electrophysiological data. 
            % Frontiers in Perception Science:233. 
            % http://dx.doi.org/10.3389/fpsyg.2012.00233.
            
            % Acunzo DJ, MacKenzie G, van Rossum MCW. 2012. 
            % Systematic biases in early ERP and ERF components as a result of high-pass filtering. 
            % Journal of Neuroscience Methods 209:212–218. 
            % http://dx.doi.org/10.1016/j.jneumeth.2012.06.011.
    
        % Practical Matlab design
        % e.g. http://www.mathworks.com/help/signal/ref/filtfilt.html        
                
        % The filtering here is done by the basicfilter() from ERPLab        
        N = filterOrder;
        lowPass = cutOffs(2);      
        highPass = cutOffs(1);        
        typef = 0; % IIF Butterworth
        chanArray = 1; % only one channel inside this function

        %% BANDPASS
        % Actually low pass and high pass with same order characteristics
        
            % call the subfunction from ERPLAB, slightly modified basicfilter()
            % modification in terms of outputs, not in regard to the
            % computation itself            
            try
                [dataOut, ferror, b, a, v, frec3dB_final, xdB_at_fx, orderx] = basicfilter_mod(data, chanArray, Fs, lowPass, highPass, N, typef);                  
            catch err
                err                
                msgboxText = ['Oops! filter is not working properly.\n Data have undefined numerical results.\n'...
                               'We strongly recommend that you change some filter parameters, for instance, decrease filter order.'];
                title = 'ERPLAB: basicfilter() error: undefined numerical results';
                error(sprintf(msgboxText), title);

                error('Have you added ERPLAB to Matlab Path?')
                
            end
            
                if length(dataOut(isnan(dataOut))) == length(dataOut)
                    warning('NaN vector returned from bandpass filter, you probably used too high order, try to reduce and what happens, auto-reduce the order by 2 now')
                    
                    N = N - 2;
                    [dataOut, ferror, b, a, v, frec3dB_final, xdB_at_fx, orderx] = basicfilter_mod(data, chanArray, Fs, lowPass, highPass, N, typef);                  
                    
                    if length(dataOut(isnan(dataOut))) == length(dataOut)
                       error('NaN vector still returned, check what is the problem!') 
                    end
                    
                end

                % Now the basic_filter constructs the transfer function, but
                % see for example: http://www.mathworks.com/matlabcentral/newsreader/view_thread/263952

                    %{
                    % Calculate the zpk values using the BUTTER function.
                    [z, p, k] = butter(N, Fc/(Fs/2));

                    % To avoid round-off errors, do not use the transfer function. Instead
                    % get the zpk representation and convert it to second-order sections.
                    [sos_var,g] = zp2sos(z, p, k);
                    Hd = dfilt.df2sos(sos_var, g);
                    %}        

                % apply additional BAND-PASS filter with more steep stop-bandpass
                % behavior
                %{
                if parameters.applySteepBandPass == 1
                    [dataOut2] = pre_bandPass_SteepDesign(dataOut, lowPass, highPass, Fs, filterOrderSteep, handles.parameters, handles);
                else
                    dataOut2 = dataOut;
                end
                %}
                dataOut2 = [];

        %% LOW and HIGH pass filters combined
        
            % this way we can use more aggressive filterOrder for the lowpass
            % filter than for the high-pass filter as empirically the N = 6 for
            % high pass filter is kinda maximum whereas you can go higher with
            % low-pass filter without getting a vector full of NaNs

                %{
                % filter_tf is again from the ERPLAB        
                [b_lo, a_lo, labelf_lo, v_lo, frec3dB_final_lo, xdB_at_fx_lo, orderx_lo] = filter_tf(typef, orderLow, lowPass, 0, Fs);
                [b_hi, a_hi, labelf_hi, v_hi, frec3dB_final_hi, xdB_at_fx_hi, orderx_hi] = filter_tf(typef, orderHigh, 0, highPass, Fs);

                % Butterworth bandpass (cascade), this also inside the
                % basicfilter_mod
                dataOut = filtfilt(b_lo, a_lo, data); % 
                dataOut = filtfilt(b_hi, a_hi, dataOut);

                plot_debugBandPassPlot(data, dataOut, dataOut2, ferror, b_lo, a_lo, b_hi, a_hi, v_lo, frec3dB_final_lo, xdB_at_fx_lo, orderx_lo, Fs, lowPass, highPass, N, typef, handles.parameters, handles)
                %}
                
        % Plot the results if wanted
        if handles.flags.showDebugPlots == 1
            plot_debugBandPassPlot(data, dataOut, dataOut2, ferror, b, a, v, frec3dB_final, xdB_at_fx, orderx, Fs, lowPass, highPass, N, typef, handles.parameters, handles)
            pause
        end
        
        % output
        y = dataOut;
       