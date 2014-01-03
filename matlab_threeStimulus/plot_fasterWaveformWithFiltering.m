function plot_fasterWaveformWithFiltering(t, averWaveForm_in, averWaveForm, noOfValidTrials_in, noOfValidTrials, erpType, parameters, handles)

    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempWaveformPlot.mat';
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
    
    style = handles.style;

    scrsz = style.scrsz;
    fig = figure('Color', 'w', 'Name', ['Waveform - ', erpType]);
        set(fig, 'Position', [0.02*scrsz(3) 0.2*scrsz(4) 0.95*scrsz(3) 0.65*scrsz(4)])
    	rows = 2;
        cols = 3;
        
    %% IN
    ind = 1;
    sp(ind) = subplot(rows,cols,ind);
    
        hold on
        p(ind,:) = plot(t, averWaveForm_in(1:parameters.EEG.nrOfChannels,:));
        pEOG = plot(t, averWaveForm_in(parameters.EEG.nrOfChannels+1:parameters.EEG.nrOfChannels+1,:), 'm');        
        lin(ind) = line([min(t) max(t)], [0 0]);
        hold off
        
            filtType = 'IN (w artifacts)';
            titStr = sprintf('%s\n%s\n%s', filtType, ['n = [', num2str(noOfValidTrials_in), ']'], [num2str(parameters.filter.bandPass_loFreq), '-', num2str(parameters.filter.bandPass_hiFreq), ' Hz (order=', num2str(parameters.filterOrder), ')']);
            tit(ind) = title(titStr);
            
            lab(ind,1) = xlabel('Time [ms]');
            lab(ind,2) = ylabel('Amplitude [\muV]');
            
            leg = legend([p, pEOG], 'Fz', 'Cz', 'Pz', 'Oz', 'EOG');
                set(leg, 'Position',[0.049 0.794 0.0512 0.128])
                legend('boxoff')
    
      
    %% OUT ('GENERAL')
    ind = ind+1;
    sp(ind) = subplot(rows,cols,ind);
    
        hold on
        p(ind,:) = plot(t, averWaveForm);
        lin(ind) = line([min(t) max(t)], [0 0]);
        hold off
        
            filtType = 'OUT (wo artifacts)';
            titStr = sprintf('%s\n%s\n%s', filtType, ['n = [', num2str(noOfValidTrials), ']'], [num2str(parameters.filter.bandPass_loFreq), '-', num2str(parameters.filter.bandPass_hiFreq), ' Hz (order=', num2str(parameters.filterOrder), ')']);
            tit(ind) = title(titStr);
            
            lab(ind,1) = xlabel('Time [ms]');
            lab(ind,2) = ylabel('Amplitude [\muV]');
    
        
    
    %% ERP
    ind = ind+1;
    sp(ind) = subplot(rows,cols,ind);   
        
        % filter
        averWaveForm_filt = zeros(size(averWaveForm,2), size(averWaveForm,1));
        loFreqCut = parameters.filter.bandPass_ERP_loFreq;
        hiFreqCut = parameters.filter.bandPass_ERP_hiFreq;
        filterOrder = parameters.filterOrder_ERP;        
        for ch = 1 : parameters.EEG.nrOfChannels
            averWaveForm_filt(:,ch) = pre_bandbassFilter(averWaveForm(ch,:)', parameters.EEG.srate, [hiFreqCut loFreqCut], filterOrder, [], handles);             
        end
    
        hold on
        p(ind,:) = plot(t, averWaveForm_filt');
        lin(ind) = line([min(t) max(t)], [0 0]);
        hold off
        
            filtType = 'ERP';
            titStr = sprintf('%s\n%s', filtType, [num2str(loFreqCut), '-', num2str(hiFreqCut), ' Hz (order=', num2str(filterOrder), ')']);
            tit(ind) = title(titStr);
            
            lab(ind,1) = xlabel('Time [ms]');
            lab(ind,2) = ylabel('Amplitude [\muV]');
    
        % define global y-limits
        yLims = [1.1*min(min(averWaveForm)) max(max(averWaveForm))*1.1];
               
            
    
    %% ALPHA
    ind = ind+1;
    sp(ind) = subplot(rows,cols,ind);
    
        % filter
        averWaveForm_filt = zeros(size(averWaveForm,2), size(averWaveForm,1));
        loFreqCut = parameters.filter.bandPass_Alpha_loFreq;
        hiFreqCut = parameters.filter.bandPass_Alpha_hiFreq;
        filterOrder = parameters.filterOrder_Alpha;        
        for ch = 1 : parameters.EEG.nrOfChannels
            averWaveForm_filt(:,ch) = pre_bandbassFilter(averWaveForm(ch,:)', parameters.EEG.srate, [hiFreqCut loFreqCut], filterOrder, [], handles); 
        end
    
        hold on
        p(ind,:) = plot(t, averWaveForm_filt');
        lin(ind) = line([min(t) max(t)], [0 0]);
        hold off
        
            filtType = 'ALPHA';
            titStr = sprintf('%s\n%s', filtType, [num2str(loFreqCut), '-', num2str(hiFreqCut), ' Hz (order=', num2str(filterOrder), ')']);
            tit(ind) = title(titStr);
            
            lab(ind,1) = xlabel('Time [ms]');
            lab(ind,2) = ylabel('Amplitude [\muV]');
    
    
    %% CNV
    ind = ind+1;
    sp(ind) = subplot(rows,cols,ind);
    
        % filter
        averWaveForm_filt = zeros(size(averWaveForm,2), size(averWaveForm,1));
        loFreqCut = parameters.filter.bandPass_CNV_loFreq;
        hiFreqCut = parameters.filter.bandPass_CNV_hiFreq;
        filterOrder = parameters.filterOrder_CNV;        
        for ch = 1 : parameters.EEG.nrOfChannels
            averWaveForm_filt(:,ch) = pre_bandbassFilter(averWaveForm(ch,:)', parameters.EEG.srate, [hiFreqCut loFreqCut], filterOrder, [], handles); 
        end
    
        hold on
        p(ind,:) = plot(t, averWaveForm_filt');
        lin(ind) = line([min(t) max(t)], [0 0]);
        hold off
        
            filtType = 'CNV';
            titStr = sprintf('%s\n%s', filtType, [num2str(loFreqCut), '-', num2str(hiFreqCut), ' Hz (order=', num2str(filterOrder), ')']);
            tit(ind) = title(titStr);
            
            lab(ind,1) = xlabel('Time [ms]');
            lab(ind,2) = ylabel('Amplitude [\muV]');
    
    
    %% P300
    ind = ind+1;
    sp(ind) = subplot(rows,cols,ind);
    
        % filter
        averWaveForm_filt = zeros(size(averWaveForm,2), size(averWaveForm,1));
        loFreqCut = parameters.filter.bandPass_P300_loFreq;
        hiFreqCut = parameters.filter.bandPass_P300_hiFreq;
        filterOrder = parameters.filterOrder_ERP;        
        for ch = 1 : parameters.EEG.nrOfChannels
            averWaveForm_filt(:,ch) = pre_bandbassFilter(averWaveForm(ch,:)', parameters.EEG.srate, [hiFreqCut loFreqCut], filterOrder, [], handles); 
        end
    
        hold on
        p(ind,:) = plot(t, averWaveForm_filt');
        lin(ind) = line([min(t) max(t)], [0 0]);
        hold off
        
            filtType = 'P300';
            titStr = sprintf('%s\n%s', filtType, [num2str(loFreqCut), '-', num2str(hiFreqCut), ' Hz (order=', num2str(filterOrder), ')']);
            tit(ind) = title(titStr);
            
            lab(ind,1) = xlabel('Time [ms]');
            lab(ind,2) = ylabel('Amplitude [\muV]');
    
            
        % annotate a note that the pre-baseline is not correct for the
        % filtered ones, however this is not a problem as these waveforms
        % here are not the final output filtered waveforms, these are just
        % for illustration
        axes(sp(cols+1))
            txAnnot = text(min(t), min(yLims) * 1.4, ['Note that the filtered waveforms are for visualization only, and are not', ...            
                        ' the final output waveforms, thus the pre-stimulus baselines are not correct as they are normalized based on the "GENERAL" ERP']);
            
    %% General styling
    set(lin, 'Color', 'k')
    
    set(sp, 'XLim', [min(t) max(t)], 'YLim', yLims)
    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold') 
    
    set(txAnnot, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'HorizontalAlignment', 'left', 'FontAngle', 'Italic')
    
    %% Auto-SAVE
    try
        if handles.figureOut.ON == 1                     
            drawnow
            dateStr = plot_getDateString(); % get current date as string          
            fileNameOut = sprintf('%s%s%s%s', 'debug_FASTERwaveform', '_', strrep(handles.inputFile, '.bdf', ''), erpType,  '.png');
            export_fig(fullfile(handles.path.debugFASTER, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            %cd(path.code)
        end
    catch err
        err
        str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
        error(str)
    end