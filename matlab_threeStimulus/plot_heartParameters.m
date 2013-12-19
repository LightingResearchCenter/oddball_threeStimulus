 function plot_heartParameters(forDebugPlot_QRS, heart, parameters, heartrateSampleRate, style, handles)
        
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempHeartPlotDebug.mat';
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
    
    scrsz = handles.style.scrsz;
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.02*scrsz(3) 0.05*scrsz(4) 0.95*scrsz(3) 0.85*scrsz(4)])
        rows = 3;
        cols = 2;
        
        zoomX = 0.4; % fraction of all x samples
        xDuration = 2.5; % seconds
        
    % ECG (QRS DETECTION)
    i = 1;
    sp(i) = subplot(rows, cols, [i i+1]);
        hold on
        p(i,1) = plot(forDebugPlot_QRS.t, forDebugPlot_QRS.ECG); 
        p(i,2:5) = plot(forDebugPlot_QRS.t_Rraw, forDebugPlot_QRS.Rraw, 'r^', ...
                        forDebugPlot_QRS.t_S, forDebugPlot_QRS.S, 'b*', ...
                        forDebugPlot_QRS.t_Q, forDebugPlot_QRS.Q, 'go', ...                        
                        forDebugPlot_QRS.t_R, forDebugPlot_QRS.R, 'k^'); 
        set(p(i,1), 'Color', [.7 .7 .7])
        hold off
        
        titStr = sprintf('%s\n%s', ['ECG (QRS Detection), f = ', num2str(heartrateSampleRate), ' Hz'], strrep(handles.inputFile, '.bdf', ''));
        tit(i) = title(titStr, 'Interpreter', 'none');
        lab(i,1) = xlabel('Time [s]'); 
        lab(i,2) = ylabel('Amplitude [mV]');        
        leg(i) = legend('ECG','R raw','S','Q', 'R');
            set(leg(i), 'Position',[0.913690476190475 0.79535274356103 0.0541979949874687 0.120380739081747])
            legend('boxoff')  
        xlim([min(forDebugPlot_QRS.t) max(forDebugPlot_QRS.t)])
        drawnow

    % RR TIME SERIES
    i = i + 1;
    sp(i) = subplot(rows, cols, i+1);    
        p(i,:) = plot(heart.vector.RR_t, heart.vector.RR_timeSeriesAfterOutlier, 'b');
        titStr = sprintf('%s\n%s', 'RR Timeseries', ['HR, mean = ', num2str(heart.scalar.HR_Mean), ', median = ', num2str(heart.scalar.HR_Median)]);
        tit(i) = title(titStr);
        lab(i,1) = xlabel('Time [s]'); 
        lab(i,2) = ylabel('RR Interval [ms]');                
        leg(i) = legend(' '); legend('hide')
        xlim([min(heart.vector.RR_t) max(heart.vector.RR_t)])
        drawnow

    % HRV SPECTRUM
    i = i + 1;
    sp(i) = subplot(rows, cols, i+1);
    
        hold on                    
        area(heart.vector.F, heart.vector.PSD_In, 'FaceColor', 'r', 'EdgeColor', 'none');
        p(i,:) = plot(heart.vector.F, heart.vector.PSD, 'b', 'LineWidth', 1);
        
        tit(i) = title('HRV Spectrum (Lomb-Scargle PSD)');
        lab(i,1) = xlabel('Frequency [Hz]');
        lab(i,2) = ylabel('Power [s^2Hz^{-1}]');
        leg(i) = legend('Filtered', 'Input');
            set(leg(i), 'Position',[0.916196741854636 0.571248600223964 0.0604636591478697 0.0510918253079507])
            legend('boxoff')  
        
        xlim([0 parameters.heart.freqBins.HF(2)])
        yLims = get(gca, 'YLim');
        
        % annotate bands
        lin(1,1) = line([parameters.heart.freqBins.LF(1) parameters.heart.freqBins.LF(1)], yLims);
        lin(1,2) = line([parameters.heart.freqBins.LF(2) parameters.heart.freqBins.LF(2)], yLims);       
        lin(2,1) = line([parameters.heart.freqBins.HF(1) parameters.heart.freqBins.HF(1)], yLims);
        lin(2,2) = line([parameters.heart.freqBins.HF(2) parameters.heart.freqBins.HF(2)], yLims);
        lin(3,1) = line([parameters.heart.freqBins.VLF(1) parameters.heart.freqBins.VLF(1)], yLims);
        lin(3,2) = line([parameters.heart.freqBins.VLF(2) parameters.heart.freqBins.VLF(2)], yLims);
            set(lin, 'Color', [.2 .2 .2], 'LineStyle', '--')
            
            offs = 0.9;
            tx(1) = text(mean(parameters.heart.freqBins.LF), yLims(2)*offs, 'LF');
            tx(2) = text(mean(parameters.heart.freqBins.HF), yLims(2)*offs, 'HF');
            tx(3) = text(mean(parameters.heart.freqBins.VLF), yLims(2)*offs, 'VLF');
            drawnow
            
        % annotate HRV parameters
        xLoc = 0.97 * parameters.heart.freqBins.HF(2);        
        yRange = yLims(2) - yLims(1);
        yInterv = yRange / 12;
        yOffs = yRange / 10;
        ind = 1; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['IBI = ', num2str(heart.scalar.avIBI), ' ms']);
        ind = 2; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['RMSSD = ', num2str(heart.scalar.RMSSD), ' ms']);
        ind = 3; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['SDNN = ', num2str(heart.scalar.SDNN), ' ms']);
        ind = 4; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['VLF = ', num2str(heart.scalar.VLF), ' s^2']);
        ind = 5; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['LF = ', num2str(heart.scalar.LF), ' s^2']);
        ind = 6; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['HF = ', num2str(heart.scalar.HF), ' s^2']);
        ind = 7; txHRV(ind) = text(xLoc, yLims(2) - yOffs - (ind*yInterv), ['LF/HF = ', num2str(heart.scalar.NLF), '']);
        
    % DFA
    i = i + 1;
    sp(i) = subplot(rows, cols, i+1);
    
        alpha = heart.scalar.DFA_alphaScaling;
        intervals = heart.vector.DFA_intervals;
        flucts = heart.vector.DFA_flucts;
        
        fittingOffsetIndex = 1; % choose of what data point has to be crossed for linear regression
        offset1 = (flucts(fittingOffsetIndex)) - (intervals(fittingOffsetIndex) .^ alpha);
        %fittingOffsetIndex = 4; % choose of what data point has to be crossed for linear regression
        alpha_2 = alpha - 1;        
        offset_2 = (flucts(fittingOffsetIndex) ./ intervals(fittingOffsetIndex)) - (intervals(fittingOffsetIndex) .^ alpha_2);
    
        pDFA(i,:) = loglog(intervals, flucts, ...
                               intervals, ((intervals .^ alpha) + offset1), 'o', ...
                               intervals, (flucts ./ intervals), ...
                               intervals, ((intervals .^ alpha_2) + offset_2), 'o');
            
        tit(i) = title('Detrended Fluctuation Analysis (DFA)');
        lab(i,1) = xlabel('Time Scale n(data points)');
        
        yLabString = sprintf('%s\n%s', 'Detrended fluctuation F(n)', 'F(n) / n');
        lab(i,2) = ylabel(yLabString);
        leg(i) = legend(['flucts, \alpha = ', num2str(alpha,2)], '((intervals ^ alpha) + offset_1)', ...
                        ['(flucts / intervals), \alpha = ', num2str(alpha_2,2)], '((intervals ^ alpha_2) + offset_2)', 'Location', 'Best');
                        legend('boxoff')  
                        
        set(pDFA(i,1:2), 'Color', [0 0.533 0.831], 'MarkerSize', 7, 'MarkerFaceColor', [0 0.533 0.831], 'MarkerEdgeColor', [0 0 0]);
        set(pDFA(i,3:4), 'MarkerSize', 7, 'MarkerFaceColor', [1 0.4 0], 'MarkerEdgeColor', [0 0 0]);
        drawnow

    % MFDFA
    i = i + 1;
    sp(i) = subplot(rows, cols, i+1);
    
        h = heart.vector.MFDFA_h;
        Dh = heart.vector.MFDFA_Dh;
    
        p = plot(h, Dh , 'ko');
        lab(i,1) = xlabel('h (Hölder/Hurst Exponent)');
        yStr = sprintf('%s\n%s', 'D(h) (FractalDimension', '/MultifractalSpectrum)');
        lab(i,2) = ylabel(yStr);
        set(p, 'markerFaceColor', 'b')
        tit(i) = title('Multifractal DFA (MFDFA)');

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
            drawnow

    
    % STYLE                
    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)                 
    set(txMFDFA, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')       
    set(txHRV, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold', 'HorizontalAlignment', 'right')
    set(tx, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')       
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold') 
    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')
    
    % Auto-SAVE
    try
        if handles.figureOut.ON == 1      
            drawnow
            dateStr = plot_getDateString(); % get current date as string          
            fileNameOut = sprintf('%s%s', 'heartSummary_', strrep(handles.inputFile, '.bdf', ''), '.png');
            export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            a = 1
            %cd(path.code)
        end
    catch err
        err
        str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
        error(str)
    end
