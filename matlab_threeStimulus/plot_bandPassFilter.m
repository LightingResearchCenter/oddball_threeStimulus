function plot_bandPassFilter(CNV, filt, raw, handles)

    debugMatFileName = 'tempBandpassFilter.mat';
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

    yLims = [-30 30];
    
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.05*scrsz(3) 0.15*scrsz(4) 0.92*scrsz(3) 0.9*scrsz(4)])  
        
        rows = 2;
        cols = 1;
        
        zoomRange = [402000 412000];
        
    i = 1;
    sp(i) = subplot(rows,cols,i);
        
        hold on
        p(i,1) = plot(CNV(:,1), 'r');
        p(i,2) = plot(filt(:,2), 'g'); 
        p(i,3) = plot(filt(:,1), 'b');
        p(i,4) = plot(CNV(:,2), 'k');
%         p(i,3) = plot(raw(:,1), 'k');
%         p(i,4) = plot(raw(:,2), 'b');
        hold off
        
            ylab(i) = ylabel('\muV');
            xlab(i) = xlabel('Sample');
            set(sp(i), 'XLim', [1 length(filt(:,2))])
            
            tit = title('Bandpass plot for CNV and ERP component analysis before EP denoising');
        
    i = 2;
    sp(i) = subplot(rows,cols,i);
        
        hold on
        p(i,1) = plot(CNV(:,1), 'r');        
        p(i,2) = plot(filt(:,2), 'g');        
        p(i,3) = plot(filt(:,1), 'b');
        p(i,4) = plot(CNV(:,2), 'k');
%         p(i,3) = plot(raw(:,1), 'k');
%         p(i,4) = plot(raw(:,2), 'b');
        hold off
        
            ylab(i) = ylabel('\muV');
            xlab(i) = xlabel('Sample');
            set(sp(i), 'XLim', zoomRange)
            
            leg = legend('for CNV', 'for ERP', '"ERP" at Cz', '"CNV" at Fz');
                legend('boxoff')
        
                yOffset = 3;
                xOffset = 100;
                tx(1) = text(zoomRange(1)+xOffset, yLims(2)-yOffset, 'Filter Characteristics');
                tx(2) = text(zoomRange(1)+xOffset, yLims(2)-yOffset*2, ['ERP Components - Bandpass: ', num2str(handles.parameters.filter.bandPass_loFreq), '-', num2str(handles.parameters.filter.bandPass_hiFreq), ' Hz (order = ', num2str(handles.parameters.filterOrder), ')']);
                tx(3) = text(zoomRange(1)+xOffset, yLims(2)-yOffset*3, ['CNV - Bandpass: ', num2str(handles.parameters.filter.bandPass_CNV_loFreq), '-', num2str(handles.parameters.filter.bandPass_CNV_hiFreq), ' Hz (order = ', num2str(handles.parameters.filterOrder_CNV), ')']);
                
                tx(4) = text(zoomRange(1)+xOffset, yLims(1)+yOffset*2, 'Red and Black should be overlapping, as well as Green and Blue, not deterministic the filter / artifact rejection?');
        
    % General styling    
    set(sp, 'YLim', yLims)
    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase)
    set([xlab ylab], 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase)
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase+2, 'FontWeight', 'bold')
    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase)
    set(tx, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase)
    set(tx(1), 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
    
    % Auto-SAVE
    try
        if handles.figureOut.ON == 1      
            drawnow
            dateStr = plot_getDateString(); % get current date as string          
            fileNameOut = sprintf('%s%s', 'bandpassFilterPlot_', strrep(handles.inputFile, '.bdf', ''), '.png');
            export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            %cd(path.code)
        end
    catch err
        err
        str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
        error(str)
    end

