function plot_artifactRegressionPlot(t, dataIn, afterEOG_removal, afterECG_removal, dataOut, ch, zoomRange, parameters, handles)

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempArtifactRegressionPlot.mat';
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

    % DEBUG difference matrices   
    EOG = dataIn(:,parameters.EEG.nrOfChannels+1);
    ECG = dataIn(:,parameters.EEG.nrOfChannels+2);

    diff_inEog = dataIn(:,1:parameters.EEG.nrOfChannels) - afterEOG_removal;
    diff_eogEcg = afterEOG_removal - afterECG_removal;
    diff_inEcg = dataIn(:,1:parameters.EEG.nrOfChannels) - afterECG_removal;
    diff_total = dataIn(:,1:parameters.EEG.nrOfChannels) - dataOut;

    % debug plot
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color','w','Name','Regression Artifacts');      
        set(fig, 'Position', [0.05*scrsz(3) 0.15*scrsz(4) 0.92*scrsz(3) 0.9*scrsz(4)])                
        rows = 4;
        cols = 2;
        
    i = 1;
    sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, diff_inEog(:,ch), t, dataIn(:,ch));
        tit(1) = title(' ');
        lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

        i = i + 1;
        sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, diff_inEog(:,ch), t, dataIn(:,ch));
        leg(1) = legend('Diff(In-EOG)', 'In');
            legend('boxoff')
            lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

    i = i + 1;
    sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, diff_inEcg(:,ch), t, diff_eogEcg(:,ch));
        tit(2) = title(' ');
        lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

        i = i + 1;
        sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, diff_inEcg(:,ch), t, diff_eogEcg(:,ch));
        leg(2) = legend('Diff(In-ECG)', 'Diff(EOG-ECG)');
            legend('boxoff')
            lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

    i = i + 1;
    sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, diff_eogEcg(:,ch), t, diff_total(:,ch));
        tit(3) = title(' ');
        lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

        i = i + 1;
        sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, diff_eogEcg(:,ch), t, diff_total(:,ch));
        leg(3) = legend('Diff(EOG-ECG)', 'Diff(In-Out)');
            legend('boxoff')
            lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

    i = i + 1;
    sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, dataIn(:,ch), t, dataOut(:,ch));
        tit(4) = title(' ');
        lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

        i = i + 1;
        sp(i) = subplot(rows,cols,i); 
        p(i,:) = plot(t, dataIn(:,ch), t, dataOut(:,ch));
        leg(4) = legend('In', 'Out');
            legend('boxoff')
            lab(i,1) = xlabel('Time [s]'); lab(i,2) = ylabel('\muV');

    % STYLE            
    set(sp, 'XLim', [min(t) max(t)])
    set(sp([2 4 6 8]), 'XLim', zoomRange)

        set(sp([1 2]), 'YLim', 60*handles.style.ERP_yLimits)
        %set(sp([9]), 'YLim', 10*handles.style.ERP_yLimits)
        %set(sp([10]), 'YLim', 5*handles.style.ERP_yLimits)

    set(leg, 'Location', 'NorthEastOutside')
    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
    
    set(p(:,1), 'Color', [0.1 0.6 1])
    set(p(:,2), 'Color', [1 0.2 0.6])

    % Auto-SAVE
    try
        if handles.figureOut.ON == 1      
            drawnow
            dateStr = plot_getDateString(); % get current date as string          
            fileNameOut = sprintf('%s%s', 'artifactRemoval_', strrep(handles.inputFile, '.bdf', ''), '.png');
            export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            %cd(path.code)
        end
    catch err
        err
        str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
        error(str)
    end