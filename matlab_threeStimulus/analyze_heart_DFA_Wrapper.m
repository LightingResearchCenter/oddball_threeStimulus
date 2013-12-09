function heart = analyze_heart_DFA_Wrapper(heart, hp, thp, rrTimes, rrPeakInterval, rrPeakAmplitude, heartrateSampleRate, parameters, handles)

    debugMatFileName = 'tempHeartDFA.mat';
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
    
    debugON = 1; 
    
    %% Default KARDIA Implementation
    
        %{
        parameters.heart.DFAparam_sliding = 's'; % adding overlap increases computational cost but gives better precision for the estimate
        parameters.heart.DFAparam_graph = 0; % plot or not inside kardia_DFA
        parameters.heart.DFAparam_minbox = 4;
        parameters.heart.DFAparam_maxbox = floor(length(hp)/4);    
        [alpha,n,Fn] = kardia_DFA(hp, parameters.heart.DFAparam_sliding, parameters.heart.DFAparam_graph, parameters.heart.DFAparam_minbox, parameters.heart.DFAparam_maxbox)

            plot(log10(n),log10(Fn))
                xlabel('log(n)')
                ylabel('log(Fn)')
                title('Detrended Fluctuation Analysis')
        %}
    
    %% FAST DFA from OXFORD
        
        disp(['              - Monofractal DFA (fastdfa)'])
        % http://www.mathworks.com/matlabcentral/fileexchange/19795-detrended-fluctuation-analysis
        % http://www.eng.ox.ac.uk/samp/dfa_soft.html 
        [alpha, intervals, flucts] = fastdfa(hp');
        
            % Outputs:
            %    alpha      - Estimated scaling exponent
            %    intervals  - List of sample interval widths at each scale
            %    flucts     - List of fluctuation amplitudes at each scale
            heart.scalar.DFA_alphaScaling = alpha;
            heart.vector.DFA_intervals = intervals;
            heart.vector.DFA_flucts = flucts;
            
    %% Multi-fractal DFA   
                           
        disp(['              - Multifractal DFA (MFDFA)'])
        
        m = 1;        
        scmin = heartrateSampleRate;        
        scmax = length(hp);
        ressc = 100;
        qmin = -5;
        qmax = 5;
        qres = abs(((qmin - qmax) / 0.2))+1;    
        % eps = 1; % 1 millisecond with HRV
        warning off 
        try
            tic                      
            [s,q,Hq,h,Dh,logFq] = MFDFA(hp,m,scmin,scmax,ressc,qmin,qmax,qres);            
            timing = toc;
        catch err     
            
            if strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                err
                error('Have you added the MFDFA with subfolders to your Matlab path, as the function was not found now!')
                
            elseif strcmp(err.identifier, 'MATLAB:badRectangle')
                
                % Improper assignment with rectangular empty matrix.
                % Error in MFDFA_eps (line 78)
                %    Fq0(ns)=exp(0.5*mean(log(RMS0{ns}.^2)));
                disp('                  -> bad rectangle ')
                err
                err.identifier       
                h = NaN;
                Dh = NaN;
            
            elseif strcmp(err.identifier, 'MATLAB:badsubscript')
                
                disp('                  -> bad susubscript ')
                h = NaN;
                Dh = NaN;
                
            else
                err
                err.identifier
                disp('                  -> Some other error here, add cases while you get new errors!')
                h = NaN;
                Dh = NaN;
            end
        end

        
        % Output
        % Variables of interest from Zorick and Mandelkern (2013)
        % http://dx.doi.org/10.1371/journal.pone.0068360
        
        % see short description of how these were calculated from the 2nd
        % column of page 2 of the Zorick and Mandelkern
        heart.scalar.MFDFA_mean_h = mean(h);
        heart.scalar.MFDFA_mean_Dh = mean(Dh);
        heart.scalar.MFDFA_width_h = max(h) - min(h);
        heart.scalar.MFDFA_height_Dh = max(Dh) - min(Dh);
        
        if debugON == 1
        
            fig = figure('Color', 'w');
            
            % as Figure 4 and Figure 5 of Zorick and Mandelkern (2013)
                % Zorick T, Mandelkern MA. 2013. 
                % Multifractal Detrended Fluctuation Analysis of Human EEG: Preliminary Investigation and Comparison with the Wavelet Transform Modulus Maxima Technique. 
                % PLoS ONE 8:e68360. 
                % http://dx.doi.org/10.1371/journal.pone.0068360
                
                p = plot(h, Dh , 'ko');
                xlabel('h (Hölder/Hurst Exponent)')
                ylabel('D(h) (FractalDimension/MultifractalSpectrum)')
                set(p, 'markerFaceColor', 'b')

                hold on
                l = line([heart.scalar.MFDFA_mean_h heart.scalar.MFDFA_mean_h], [0 1]);
                hold off
                
                xlims = get(gca, 'XLim');
                ylims = get(gca, 'YLim');
                yOffset = 0.1;
                
                tx(1) = text(0.95*xlims(2), ylims(2)*(1-yOffset*1), ['mean_h = ', num2str(heart.scalar.MFDFA_mean_h)]);
                tx(2) = text(0.95*xlims(2), ylims(2)*(1-yOffset*2), ['mean_D(h) = ', num2str(heart.scalar.MFDFA_mean_Dh)]);
                tx(3) = text(0.95*xlims(2), ylims(2)*(1-yOffset*3), ['width_h = ', num2str(heart.scalar.MFDFA_width_h)]);
                tx(4) = text(0.95*xlims(2), ylims(2)*(1-yOffset*4), ['height_D(h) = ', num2str(heart.scalar.MFDFA_height_Dh)]);
                set(tx, 'HorizontalAlignment', 'right')
                
        end        
       
            
            
        % Matlab code 14------------------------------------------
        %{
        Ht_row=Ht(:);
        BinNumb=round(sqrt(length(Ht_row)));
        [freq,Htbin]=hist(Ht_row,BinNumb);
        Ph=freq./sum(freq);
        Ph_norm=Ph./max(Ph);
        Dh=1-(log(Ph_norm)./-log(mean(scale)));
        plot13;
        %}
        
        % Multifractal detrended fluctuation analysis (MFDFA). This Matlab function
        % estimate the multifractal spectrum Dh directly without q-order statistics.
        %
        % [Ht,Htbin,Ph,Dh] = MFDFA2(signal,scale,m,Fig);
        %
        % INPUT PARAMETERS---------------------------------------------------------
        %
        % signal:       input signal
        % scale:        vector of scales 
        % m:            polynomial order for the detrending
        % Fig:          1/0 flag for output plot of F0, Ht, Ph, and Dh.
        %
        % OUTPUT VARIABLES---------------------------------------------------------
        %
        % Ht:           Time evolution of the local Hurst exponent
        % Htbin:        Bin senters for the histogram based estimation of Ph and Dh 
        % Ph:           Probability distribution of the local Hurst exponent Ht
        %               is the same as Hölder exponent
        % Dh:           Multifractal spectrum 
        %
        % EXAMPLE------------------------------------------------------------------
        %
        % load fractaldata
        % scale=[7,9,11,13,15,17];
        % m=2;
        % signal1=multifractal;
        % signal2=monofractal;
        % signal3=whitenoise;
        % [Ht1,Htbin1,Ph1,Dh1] = MFDFA2(signal1,scale,m,1);
        % [Ht2,Htbin2,Ph2,Dh2] = MFDFA2(signal2,scale,m,1);
        % [Ht3,Htbin3,Ph3,Dh3] = MFDFA2(signal3,scale,m,1);
        %--------------------------------------------------------------------------
