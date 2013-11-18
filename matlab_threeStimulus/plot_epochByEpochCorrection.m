function plot_epochByEpochCorrection(ij, fig, epochIn, epochCorr, epochMean, dataType, condition, parameters, handles)        
                
    [rows,cols] = size(epochIn);
 
    t = 1000 * (((linspace(1,rows,rows))' / handles.parameters.EEG.srate) - parameters.oddballTask.ERP_baseline);

    baseline = 0;
    yOffset = 20; % [uV]
    chOffset = 2;        

    spRows = 6;
    spCols = 1;

    % EEG & EOG
    sp(1) = subplot(spRows, spCols, [1:5]);
    for i = 1 : cols-1
        p(i,1) = plot(t,epochIn(:,i) + (i*yOffset), 'b');
        hold on
        p(i,2) = plot(t,epochCorr(:,i) + (i*yOffset), 'k');            
        p(i,3) = plot(t,epochMean + (i*yOffset), 'r');
        p(i,4) = line([min(t) max(t)], [i*yOffset i*yOffset], 'Color', [.3 .3 .3]);
        yTickPos(i) = i*yOffset;
        yTickLabel{i} =  parameters.BioSemi.chName{i+2};
    end
    hold off

    titleStr = sprintf('%s\n%s\n%s', 'Baseline Removal before ICA (runica)',...
                                     'EEGLAB function: [corr,mean] = rmbase(data)',...
                                     'Note! Oz and Pz not connected to head, so they are pure noise');
        tit = title(titleStr);        

    leg = legend([p(1,1) p(1,2) p(1,3)], 'Input', 'Baseline corr', 'Mean');
        legend('boxoff')

    % ECG
    sp(2) = subplot(spRows, spCols, 6);
    i = i +1;                    
        p(i,1) = plot(t,epochIn(:,i), 'b');            
        hold on
        p(i,2) = plot(t,epochCorr(:,i), 'k');                        
        p(i,3) = plot(t,epochMean, 'r');
        p(i,4) = line([min(t) max(t)], [0 0], 'Color', [.3 .3 .3]);
        hold off
        yTickPos2 = 0;
        yTickLabel2 = parameters.BioSemi.chName{i+2};

    xlab = xlabel([condition, ', ', dataType, ', epoch #', num2str(ij)]);

    % style
    set(sp(1), 'YTick', yTickPos, 'YTickLabel', yTickLabel)    
    set(sp(2), 'YTick', yTickPos2, 'YTickLabel', yTickLabel2)    
    set(p(:,2), 'LineWidth', 2)
    set(sp, 'XLim', [min(t) max(t)])
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase+2, 'FontWeight', 'bold')
    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase)
    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase)
    set(xlab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase+1, 'FontWeight', 'bold')

    % Auto-SAVE
    try
        if handles.figureOut.ON == 1  && handles.flags.showDebugPlots == 1
            drawnow
            dateStr = plot_getDateString(); % get current date as string
            %cd(path.outputFigures)   

            % correct i string
            if i < 10
                iStr = sprintf('%s%s', num2str(0), num2str(ij));
            else
                iStr = num2str(ij);
            end

            fileNameOut = sprintf('%s%s%s%s', 'plot_epochByEpoch_', strrep(handles.inputFile, '.bdf', ''), ...
                '_', dataType, '_', condition, '_epoch',  iStr, '.png');
            export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
            %cd(path.code)
        end
    catch err
        err
        str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                      'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
        error(str)
    end

        
    