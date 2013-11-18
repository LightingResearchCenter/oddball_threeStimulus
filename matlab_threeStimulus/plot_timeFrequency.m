function plot_timeFrequency(analyzed_TF, handles)

    %% DEBUG
    debugMatFileName = 'tempTF_Plot.mat';
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
    
    chList = [1 2 3 4]; % all channels
    parameters = handles.parameters;
    style = handles.style;
    
    %% PLOT    
    
        scrsz = handles.style.scrsz;
        fig = figure('Color', 'w');
            set(fig, 'Position', [0.1*scrsz(3) 0.2*scrsz(4) 0.85*scrsz(3) 0.65*scrsz(4)])

            rows = length(chList);
            cols = 3;


            for j = 1 : length(chList)

                ch = chList(j);
                [s(j,:,:), p(j,:,:), yLims(j,:,:)] = plot_waveformOneChannel(analyzed_TF, j, rows, cols, ch, length(chList), parameters, style);
                drawnow

            end
        
    
    %% Auto-SAVE
    
        try
            
            if handles.figureOut.ON == 1      
                
                drawnow
                dateStr = plot_getDateString(); % get current date as string        
                
                loString = num2str(parameters.filter.bandPass_ERP_loFreq);
                
                if parameters.filter.bandPass_ERP_hiFreq < 10
                    hiString = ['0', num2str(parameters.filter.bandPass_ERP_hiFreq)];
                else
                    hiString = num2str(parameters.filter.bandPass_ERP_hiFreq);
                end
                
                fileNameOut = sprintf('%s%s%s%s %s%s%s%s %s%s%s%s %s', 'timeFreq_', ...                    
                    '_f0Lo-', loString, 'Hz', ...
                    '_f0Hi-', hiString, 'Hz', ...                    
                    '_EP-Scale-', num2str(parameters.ep_den.scales_postStim), ... 
                    '_', strrep(handles.inputFile, '.bdf', ''), '.png');
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
                
            end
            
        catch err
            err
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        end

    
    function [s, p, yLims] = plot_waveformOneChannel(analyzed_TF, j, rows, cols, ch, noOfChs, parameters, style)              
        
        % real part
        TF_mat_target.real.mean = analyzed_TF.target.real.mean(:,:,ch);
        TF_mat_target.real.SD = analyzed_TF.target.real.SD(:,:,ch);
        TF_mat_distracter.real.mean = analyzed_TF.distr.real.mean(:,:,ch);
        TF_mat_distracter.real.SD = analyzed_TF.distr.real.SD(:,:,ch);
        TF_mat_standard.real.mean = analyzed_TF.std.real.mean(:,:,ch);
        TF_mat_standard.real.SD = analyzed_TF.std.real.SD(:,:,ch);

        zValues = [TF_mat_target.real.mean TF_mat_distracter.real.mean TF_mat_standard.real.mean];
        yLims = [min(min(zValues)) max(max(zValues))];
        
        % imaginary part
        timeIn = analyzed_TF.distr.timep;
        freqIn = analyzed_TF.distr.freq;        
        [TIME,FREQ] = meshgrid(timeIn, freqIn);
        
        i = 1;
        ind = (j-1)*cols+i;
            s(i) = subplot(rows, cols, ind);
            p.(['ind', num2str(i)])(:,:) = plot_timeFreq_lowLevel(TIME, FREQ, TF_mat_target, parameters);

            if j == 1
                titStr = sprintf('%s\n%s\n%s', ['TARGET (', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')'], ...
                                ['freq: ', num2str(parameters.oddballTask.targetFreq), ' Hz'], ...
                                ['n = ', num2str(analyzed_TF.target.n), ' / ', num2str(parameters.oddballTask.repeatsOfCycle)]);
                tit(i) = title(titStr);
                
            else
                tit(i) = title(['(', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')']);
            end

            if j == noOfChs
                xLab(i) = xlabel('Time [ms]');
            else
                xLab(i) = xlabel('');
            end            
            yLab(i) = ylabel('Hz');
            
            %set(gca, 'XTick', timeIn)
            %set(gca, 'YTick', freqIn)

        i = 2;
        ind = (j-1)*cols+i;
            s(i) = subplot(rows, cols, ind);
            hold on
            p.(['ind', num2str(i)])(:,:) = plot_timeFreq_lowLevel(TIME, FREQ, TF_mat_distracter, parameters);

            if j == 1
                titStr = sprintf('%s\n%s\n%s', ['DISTRACTER (', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')'], ...
                                ['freq: ', num2str(parameters.oddballTask.distractFreq), ' Hz'], ...
                                [' n = ', num2str(analyzed_TF.distr.n), ' / ', num2str(parameters.oddballTask.repeatsOfCycle)]);
                tit(i) = title(titStr);
            else
                tit(i) = title(['(', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')']);
            end

            if j == noOfChs
                xLab(i) = xlabel('Time [ms]');
            else
                xLab(i) = xlabel('');
            end                
            % yLab(i) = ylabel('\muV');

        i = 3;
        ind = (j-1)*cols+i;
            s(i) = subplot(rows, cols, ind);
            hold on
            p.(['ind', num2str(i)])(:,:) = plot_timeFreq_lowLevel(TIME, FREQ, TF_mat_standard, parameters);

            if j == 1
                titStr = sprintf('%s\n%s\n%s', ['STANDARD (', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')'], ...
                                ['freq: ', num2str(parameters.oddballTask.standardFreq), ' Hz'], ...
                                [' n = ', num2str(analyzed_TF.std.n), ' / ', num2str(parameters.oddballTask.nrOfStandardsPerCycle*parameters.oddballTask.repeatsOfCycle)]);
                tit(i) = title(titStr);
            else
                tit(i) = title(['(', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')']);
            end

            if j == noOfChs
                xLab(i) = xlabel('Time [ms]');
            else
                xLab(i) = xlabel('');
            end
            % yLab(i) = ylabel('\muV');            

            
    function p = plot_timeFreq_lowLevel(TIME, FREQ, TF_mat, parameters, handles)
        
        
        handles.parameters.plot.timeFreq_contourLevels = 64;
        
        meaaan = TF_mat.real.mean
        
        whos
        
        p = contourf(TIME, FREQ, TF_mat.real.mean, 64, 'EdgeColor', 'none'); 
        colorbar

        
        