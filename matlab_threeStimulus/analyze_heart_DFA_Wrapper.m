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
        [alpha, intervals, flucts] = fastdfa(hp);
        
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
        qres = abs(((qmin - qmax) / 0.2))+1; % e.g. 51
        % eps = 1; % 1 millisecond with HRV        
        try
            tic                      
            [s,q,Hq,h,Dh,logFq] = MFDFA(hp,m,scmin,scmax,ressc,qmin,qmax,qres)
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
                
                disp('                  -> bad subscript (no MFDFA returned for RR interval vector (why actually)')
                
                    % Only one value in variable "s", others are NaN, that
                    % is why the polynomial fit fails as number of data
                    % points is 1 which is the same as polynomial degree
                    
                    % Warning: Polynomial is not unique; degree >= number of data points. 
                    % > In polyfit at 71
                    % In MFDFA at 159
                    % In analyze_heart_DFA_Wrapper at 61 
                    
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
        
            scrsz = handles.style.scrsz;
            fig = figure('Color', 'w', 'Name', 'Fractal ECG');
                set(fig, 'Position', [0.5*scrsz(3) 0.16*scrsz(4) 0.47*scrsz(3) 0.77*scrsz(4)])
                rows = 3;
                cols = 1;

                ind = 1;
                sp(ind) = subplot(rows,cols,ind);
                
                    offset1 = 0;
                    pDFA(1,:) = loglog(intervals, flucts, intervals, ((intervals .^ alpha) + offset1), 'o');
                    tit(ind) = title('DFA');
                    lab(ind,1) = xlabel('Time Scale n(data points)');
                    lab(ind,2) = ylabel('Detrended fluctuation F(n)');
                    leg(ind) = legend('flucts', '(intervals .^ alpha)', 'Location', 'Best');
                        legend('boxoff')                    
                    
                    % annotate alpha
                    tx(ind) = text(max(intervals)/100, max(flucts)/2, ['\alpha = ', num2str(alpha,2)]);

                
                ind = ind+1;
                sp(ind) = subplot(rows,cols,ind);
                
                    alpha_2 = alpha - 1;
                    fittingOffsetIndex = 4; % choose of what data point has to be crossed for linear regression
                    offset_2 = (flucts(fittingOffsetIndex) ./ intervals(fittingOffsetIndex)) - (intervals(fittingOffsetIndex) .^ alpha_2);
                
                    pDFA(2,:) = loglog(intervals, (flucts ./ intervals), intervals, ((intervals .^ alpha_2) + offset_2), 'o');
                    lab(ind,1) = xlabel('Time Scale n(data points)');
                    lab(ind,2) = ylabel('F(n) / n');
                    tit(ind) = title(' ');
                    leg(ind) = legend('(flucts ./ intervals)', '((intervals .^ alpha_2) + offset_2)', 'Location', 'Best');
                        legend('boxoff')
                    
                    % annotate alpha
                    tx(ind) = text(max(intervals)/100, max(flucts ./ intervals)/2, ['\alpha" = ', num2str(alpha_2,1)]);
                    
                ind = ind+1;
                sp(ind) = subplot(rows,cols,ind);
                
                    % as Figure 4 and Figure 5 of Zorick and Mandelkern (2013)
                    % Zorick T, Mandelkern MA. 2013. 
                    % Multifractal Detrended Fluctuation Analysis of Human EEG: Preliminary Investigation and Comparison with the Wavelet Transform Modulus Maxima Technique. 
                    % PLoS ONE 8:e68360. 
                    % http://dx.doi.org/10.1371/journal.pone.0068360

                    p = plot(h, Dh , 'ko');
                    lab(ind,1) = xlabel('h (Hölder/Hurst Exponent)');
                    yStr = sprintf('%s\n%s', 'D(h) (FractalDimension', '/MultifractalSpectrum)');
                    lab(ind,2) = ylabel(yStr);
                    set(p, 'markerFaceColor', 'b')
                    tit(ind) = title('MFDFA');

                    hold on
                    l = line([heart.scalar.MFDFA_mean_h heart.scalar.MFDFA_mean_h], [0 1]);
                    hold off

                    xlims = get(gca, 'XLim');
                    ylims = get(gca, 'YLim');
                    yOffset = 0.1;

                    txMFDFA(1) = text(0.95*xlims(2), ylims(2)*(1-yOffset*1), ['mean_h = ', num2str(heart.scalar.MFDFA_mean_h)]);
                    txMFDFA(2) = text(0.95*xlims(2), ylims(2)*(1-yOffset*2), ['mean_D(h) = ', num2str(heart.scalar.MFDFA_mean_Dh)]);
                    txMFDFA(3) = text(0.95*xlims(2), ylims(2)*(1-yOffset*3), ['width_h = ', num2str(heart.scalar.MFDFA_width_h)]);
                    txMFDFA(4) = text(0.95*xlims(2), ylims(2)*(1-yOffset*4), ['height_D(h) = ', num2str(heart.scalar.MFDFA_height_Dh)]);
                    set(txMFDFA, 'HorizontalAlignment', 'right')
                    
                % STYLE
                set(pDFA(:,1), 'Color', [0 0.533 0.831], 'LineWidth', 1.5);
                set(pDFA(:,2), 'MarkerSize', 7, 'MarkerFaceColor', [1 0.4 0], 'MarkerEdgeColor', [0 0 0]);
                
                % STYLE                
                set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)                 
                set(txMFDFA, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')       
                set(tx, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')       
                set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold') 
                set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
                set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')      
                
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
