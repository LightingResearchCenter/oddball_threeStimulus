function batch_plot_subjectERPs_withAverage(conditions, sessions, erpResponses, erpTypes, filterTypes, statFields, statParam, ...
                                       noOfDataSamplesPerERP, noOfChannels, noOfSubjects, ...
                                       t, epochsMatrix, parameters, handles)

    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    
    fig = figure('Color', 'w', 'Name', ['ERP Waveforms, ', conditions{1}]);
        set(fig, 'Position', [0.05*scrsz(3) 0.20*scrsz(4) 0.92*scrsz(3) 0.62*scrsz(4)])
    
        % Customize the color scale, use the distinguishable_colors from Matlab
        % FileExchange by Tim Holy, http://www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
        % See also: http://blogs.mathworks.com/pick/2008/08/15/colors-for-your-multi-line-plots/    
        ColorSet = distinguishable_colors(noOfSubjects);
        set(gca, 'ColorOrder', ColorSet);

        rows = 3; % erp responses, target / distracter / standard
        cols = 5; % channels + mean
        
        accum = 0;
        
        for resp = 1 : length(erpResponses)
            
            erpResponseMatrix = epochsMatrix.(conditions{1}).(sessions{1}).(erpResponses{resp}).(erpTypes{1}).(filterTypes{1});
            
            statMatrix = nanmean(erpResponseMatrix.(statParam)(:,:,:),3);
            
            for ch = 1 : noOfChannels
                
                ind = (resp-1)*cols + ch;
                
                sp(ind) = subplot(rows,cols,ind);
                plot(t, squeeze(erpResponseMatrix.(statParam)(:,ch,:)))
                tit(ind) = title(parameters.BioSemi.chName{ch+parameters.BioSemi.chOffset});
                
                if ch == 1
                   accum = accum + 1;
                   yStr = sprintf('%s\n%s', erpResponses{resp}, 'Amplitude [\muV]');
                   labY(accum) = ylabel(yStr);
                end
                
                if resp == length(erpResponses)
                    labX(ch) = xlabel('Time [ms]');
                end
                
                % MEAN
                if ch == noOfChannels
                    sp(ind+1) = subplot(rows,cols,ind+1);
                    plot(t, statMatrix)
                    leg(resp) = legend('Fz', 'Cz', 'Pz', 'Oz');
                        xOffset = 0.14;
                        gcaPos = get(gca, 'Position');
                        set(leg(resp), 'Position', [gcaPos(1)+xOffset gcaPos(2) 0.0330 0.090])
                    tit(ind+1) = title('Average ERPs');
                    ylim([-20 20])
                end
                
            end           
            
            
        end        
        
        % STYLE
        set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
        set(sp, 'XLim', [min(t) max(t)])
        set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
        set(labY, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold')            
        set(labX, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold')       
            
        % AUTO-SAVE FIGURE
        try
            if handles.figureOut.debugON == 1 
                drawnow
                dateStr = plot_getDateString(); % get current date as string
                %cd(path.outputFigures)            
                fileNameOut = ['waveformERP_allChannels_', dateStr];

                disp(['         ... saving figure to disk (', fileNameOut, '.png]'])
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig(1))
                %cd(path.code)
            end
        catch err
            err
        end
    