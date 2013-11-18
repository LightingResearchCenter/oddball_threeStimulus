function plot_indivSubject(timeVector, epochs_ep, epochs, contourMode, parameters, style, handles)

    debugMatFileName = 'tempIndivSubjPlot.mat';
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

    % SUBPLOT LAYOUT
    rows = 5;
    cols = parameters.oddballTask.numberOfIrrTrialsPerCycle;
    scrsz = style.scrsz;
    
    ch = 3;

    %{
    % PLOT - version 1
    
        % IRREGULAR
        fig(1) = figure('Color', 'w');
            p{1} = plot_routineForSingleTrialERPs_1(timeVector, ch, epochs_ep.oddball_irregular, epochs.oddball_irregular, parameters, plotParam, fig(1), rows, cols, handles);    

        % REGULAR
        fig(2) = figure('Color', 'w');
            p{2} = plot_routineForSingleTrialERPs_1(timeVector, ch, epochs_ep.oddball_regular, epochs.oddball_regular, parameters, plotParam, fig(2), rows, cols, handles);
    %}        
    
    % PLOT - version 2
    
        fig2 = figure('Color', 'w');
            set(fig2, 'Position', [0.3*scrsz(3) 0.1*scrsz(4) 0.66*scrsz(3) 0.85*scrsz(4)])
            p2{1} = plot_routineForSingleTrialERPs_2(timeVector, ch, epochs_ep.oddball_regular, epochs.oddball_regular, ... 
                                                    epochs_ep.oddball_irregular, epochs.oddball_irregular, ...
                                                    contourMode, parameters, style, fig2, rows, cols, handles);    
        

    
        
    function p = plot_routineForSingleTrialERPs_2(timeVector, ch, epochsFilt_reg, epochsRaw_reg, epochsFilt_irreg, epochsRaw_irreg, contourMode, parameters, style, fig, rows, cols, handles)
        
        % Subplot layout
        rows = 4;
        cols = 2;
        
        t = timeVector * 1000;   
        n = 32; % number of trials
        chName = parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset};
        
        [noTrials, regMatrix, irregMatrix, averRaw_reg, averFilt_reg, averRaw_irreg, averFilt_irreg] = ...
            plot_ERP_constructAverERPsAndContour(epochsRaw_reg, epochsRaw_irreg, epochsFilt_reg, epochsFilt_irreg, ch, contourMode, parameters, handles);
                
        % Plot
        i = 1;
        s(i) = subplot(rows,cols,[1 3 5]);            
            colorBarLimits = [];
            [c,h(1),tit(1),xlab(i),ylab(i),zlab(i)] = plot_ERP_asContours(t, noTrials, regMatrix, n, chName, contourMode, colorBarLimits, style, parameters, handles);
        
        i = 2;
        s(i) = subplot(rows,cols,[2 4 6]);                   
            colorBarLimits = [];
            [c,h(2),tit(2),xlab(i),ylab(i),zlab(i)] = plot_ERP_asContours(t, noTrials, irregMatrix, n, chName, contourMode, colorBarLimits, style, parameters, handles);
        
        i = 3;        
        s(i) = subplot(rows,cols,7);        
            [p(i, 1:2), xlab(i), ylab(i), leg(1)] = plot_ERP_averageWaveform(t, averRaw_reg, averFilt_reg);                      
        
        i = 4;        
        s(i) = subplot(rows,cols,8);        
            [p(i, 1:2), xlab(i), ylab(i), leg(2)] = plot_ERP_averageWaveform(t, averRaw_irreg, averFilt_irreg);
            
        % Style 
        set(h,'Edgecolor','none')
        set(s, 'XLim', [min(t) max(t)])        
        set(s, 'FontName', style.fontName, 'FontSize', style.fontSizeBase)
        
        spacing = 1;
        yTicks = 1:spacing:n;
        set(s(1:2), 'YTick', yTicks)
        set(s(3:4), 'YLim', [-6 6])
        
        set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase+2, 'FontWeight', 'bold')
        set(ylab, 'FontName', style.fontName, 'FontSize', style.fontSizeBase, 'FontWeight', 'bold')
        set(xlab, 'FontName', style.fontName, 'FontSize', style.fontSizeBase, 'FontWeight', 'bold')
        set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-2, 'Location', 'NorthEast')
        
        set(p(3:4, 1), 'Color', [0.871 0.490 0]); % raw
        set(p(3:4, 2), 'Color', [0.122 0.463 1]); % filtered/denoised
       
        % Auto-SAVE
        try
            if handles.figureOut.ON == 1      
                drawnow
                dateStr = plot_getDateString(); % get current date as string          
                fileNameOut = sprintf('%s%s', 'plot_indivSubject_', contourMode, '_', strrep(handles.inputFile, '.bdf', ''), '.png');
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
            end
        catch err
            err
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        end

    %{    
    function p = plot_routineForSingleTrialERPs_1(timeVector, ch, epochsFilt, epochsRaw, parameters, style, fig, rows, cols, handles)
        
        t = timeVector;
        
        averRaw = zeros(length(epochsRaw{1}),1);
        averFilt = zeros(length(epochsFilt{1}),1);
        
        % Single-trial ERPs
        for i = 1 : length(epochsFilt)                    
            averRaw = averRaw + epochsRaw{i}(:,ch);
            averFilt = averFilt + epochsFilt{i}(:,ch);            
        end
                
        % Average ERP
        averRaw = averRaw / i;
        averFilt = averFilt / i;
        i = i + 1;   
        
            s(i) = subplot(rows,cols,[i:cols*rows]);
                p(i, 1:2) =  plot(t, averRaw, 'k', t, averFilt, 'b');
                    tit(i) = title(['Average ERP']);
                    legend('Raw', 'Denoised')
                        legend('boxoff')
                    
        
        % Style
        set(s, 'XLim', [0 max(t)])
    %}