function batch_plot_ERPs_perSomeParameter(rowParameter, colParameter, gcaParameter, ...
                                        rowVaritions, sessions, erpResponses, erpTypes, filterTypes, statFields, statParam, ...
                                        noOfDataSamplesPerERP, noOfChannels, noOfSubjects, ...
                                        t, epochsMatrix, parameters, handles)
      
    scrsz = get(0,'ScreenSize'); % get screen size for plotting
                                    
    fig = figure('Color', 'w', 'Name', ['ERP Average Waveform Comparison']);
        set(fig, 'Position', [0.075*scrsz(3) 0.15*scrsz(4) 0.92*scrsz(3) 0.65*scrsz(4)])        
    
        % Customize the color scale, use the distinguishable_colors from Matlab
        % FileExchange by Tim Holy, http://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
        % See also: http://blogs.mathworks.com/pick/2008/08/15/colors-for-your-multi-line-plots/    
        ColorSet = distinguishable_colors(noOfSubjects);
        set(gca, 'ColorOrder', ColorSet);

        rows = 3; % 
        cols = 4; % channels + mean
        
        for rowVar = 1 :  length(rowParameter)
        
            for colVar = 1 : noOfChannels

                ind = (rowVar - 1) * noOfChannels + colVar;
                sp(rowVar,colVar) = subplot(rows, cols, ind);

                    hold on
                    for resp = 1 : length(erpResponses)
                        
                        statMatrix = nanmean(epochsMatrix.(rowVaritions{rowVar}).(sessions{1}).(erpResponses{resp}).(erpTypes{1}).(filterTypes{1}).(statParam),3);
                        p(rowVar, colVar, resp) = plot(t, statMatrix(:,colVar));
                    
                    end
                    hold off
                    
                    labX(rowVar,colVar) = xlabel('Time [ms]');
                    labY(rowVar,colVar) = ylabel('Amplitude [\muV]');
                    
                    chPartStr = parameters.BioSemi.chName{colVar+parameters.BioSemi.chOffset};
                    if rowVar == 1
                        if colVar == 1
                            titString = sprintf('%s\n%s', rowVaritions{rowVar}, chPartStr);
                        elseif colVar == 2
                            titString = sprintf('%s\n%s', sessions{1}, chPartStr);
                        elseif colVar == 3
                            titString = sprintf('%s\n%s', erpTypes{1}, chPartStr);
                        else
                            titString = sprintf('%s\n%s', filterTypes{1}, chPartStr);
                        end
                    else
                        if colVar == 1
                            titString = sprintf('%s\n%s', rowVaritions{rowVar}, chPartStr);
                        else
                            titString = sprintf('%s', chPartStr);
                        end
                    end
                    tit(rowVar,colVar) = title(titString);
                    
                    if colVar == noOfChannels
                        leg(resp) = legend('Target', 'Distracter', 'Standard');
                        xOffset = 0.18;
                        gcaPos = get(gca, 'Position');
                        set(leg(resp), 'Position', [gcaPos(1)+xOffset gcaPos(2) 0.0330 0.090])
                        legend('boxoff')
                    end
            end

        end
        
        set(sp, 'YLim', [-20 20])
        set(sp, 'XLim', [min(t) max(t)])
        
        set(p(:, :, 1), 'Color', [1 0.2 0.6])
        set(p(:, :, 2), 'Color', [0 0.8 1])
        set(p(:, :, 3), 'Color', [0 0 0])
        set(p, 'LineWidth', 2)
        
        set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)         
        set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
        set(labY, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold')            
        set(labX, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold')       
        
        % AUTO-SAVE FIGURE
        try
            if handles.figureOut.debugON == 1 
                drawnow
                dateStr = plot_getDateString(); % get current date as string
                %cd(path.outputFigures)            
                fileNameOut = ['waveformERP_subjectAverages_', dateStr];

                disp(['         ... saving figure to disk (', fileNameOut, '.png]'])
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig(1))
                %cd(path.code)
            end
        catch err
            err
        end
        %}