function plot_waveforms(epochs_target, epochs_distracter, epochs_standard, ...
             epochsEP_target, epochsEP_distracter, epochsEP_standard, ...
             analyzed_target, analyzed_distracter, analyzed_standard, handles)
         
     %% DEBUG
    debugMatFileName = 'tempTrendPlot.mat';
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

    %%     
    chList = [1 2 3 4]; % all channels
    parameters = handles.parameters;
    style = handles.style;
    
    %{
    epochsEP_target
    epochsEP_target.ERP
    epochsEP_target.ERP.ERP
    epochsEP_target.ERP.ERP{1}
    epochsEP_target.ERP.ERP{1}(:,ch)
    %}
 
    
    %% CONDITION THE DATA
    
        for j = 1 : length(chList)
            
            ch = chList(j);
    
            %% CLASSICAL FILTERING (artifact removal, and bandpass filtering)
            epochs_target_ERPmat{j} = plot_epochsFromCellToMatrix(epochs_target.ERP.ERP, ch, handles);
            epochs_distracter_ERPmat{j} = plot_epochsFromCellToMatrix(epochs_distracter.ERP.ERP, ch, handles);
            epochs_standard_ERPmat{j} = plot_epochsFromCellToMatrix(epochs_standard.ERP.ERP, ch, handles);

            targetERP{j} = plot_calculateStatsOfEpochs(epochs_target_ERPmat{j}, ch, handles);
            distracterERP{j} = plot_calculateStatsOfEpochs(epochs_distracter_ERPmat{j}, ch, handles);
            standardERP{j} = plot_calculateStatsOfEpochs(epochs_standard_ERPmat{j}, ch, handles);

            %% EP DENOISED
            epochsEP_target_ERPmat{j} = plot_epochsFromCellToMatrix(epochsEP_target.ERP.ERP, ch, handles);
            epochsEP_distracter_ERPmat{j} = plot_epochsFromCellToMatrix(epochsEP_distracter.ERP.ERP, ch, handles);
            epochsEP_standard_ERPmat{j} = plot_epochsFromCellToMatrix(epochsEP_standard.ERP.ERP, ch, handles);

            targetERP_EP{j} = plot_calculateStatsOfEpochs(epochsEP_target_ERPmat{j}, ch, handles);
            distracterERP_EP{j} = plot_calculateStatsOfEpochs(epochsEP_distracter_ERPmat{j}, ch, handles);
            standardERP_EP{j} = plot_calculateStatsOfEpochs(epochsEP_standard_ERPmat{j}, ch, handles);
            
        end
    
    %% PLOT
    
    scrsz = handles.style.scrsz;
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.1*scrsz(3) 0.2*scrsz(4) 0.85*scrsz(3) 0.65*scrsz(4)])
        rows = 4;
        cols = 3;
    
        x = linspace(-length(targetERP_EP{j}.mean)/2, length(targetERP_EP{j}.mean)/2, length(targetERP_EP{j}.mean))';
        x = 1000 * x / handles.parameters.EEG.srate;
        
        alpha = 0.3; % transparency of the std shading (alpha)
        
        for j = 1 : length(chList)
        
            ch = chList(j);
            [s(j,:,:), p(j,:,:)] = plot_waveformOneChannel(x, targetERP{j}, targetERP_EP{j}, distracterERP{j}, distracterERP_EP{j}, standardERP{j}, standardERP_EP{j}, j, rows, cols, alpha, ch, length(chList), parameters, style);
        
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
                
                fileNameOut = sprintf('%s%s%s%s %s%s%s%s %s%s%s%s %s', 'waveform_', ...                    
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

            
    %% SUBFUNCTIONS
    
        function matrixOut = plot_epochsFromCellToMatrix(cellIn, ch, handles)

            numberOfEpochs = length(cellIn);
            numberOfSamples = length(cellIn{1}(:,ch));

            matrixOut = zeros(numberOfSamples, numberOfEpochs);
                % so as meny vectors as there are epochs, e.g. 40 for targets

            for i = 1 : numberOfEpochs            
                matrixOut(:,i) = cellIn{i}(:,ch);           
            end

        function statsOut = plot_calculateStatsOfEpochs(matrixIn, ch, handles)

            dimension = 2; % 40 columns -> 1 column, rows as many as there are samples
            flag = 0; % flag is 0 or 1 to specify normalization by n – 1 or n, respectively, 
                      % where n  is the number of remaining observations after removing observations with NaN values.

            NanCount = sum(isnan(matrixIn),1); % return number of NaNs per column            
            NanCountBoolean = logical(NanCount);
            
            % if a column (epoch) contains even one NaN, it is an
            % artifacted one so we need to transform the whole column to
            % NaNs (this could have been done previously also with the
            % epoching part)
            matrixIn(:,NanCountBoolean) = NaN;
                                  
            statsOut.mean = nanmean(matrixIn, dimension);
            statsOut.std = nanstd(matrixIn, flag, dimension);
            
            % n is the number of non-NaN epochs
            statsOut.n = sum(NanCountBoolean == 0);            
            
            
        function [s,p] = plot_waveformOneChannel(x, targetERP, targetERP_EP, distracterERP, distracterERP_EP, standardERP, standardERP_EP, j, rows, cols, alpha, ch, noOfChs, parameters, style)

            i = 1;
            ind = (j-1)*cols+i;
                s(i) = subplot(rows, cols, ind);
                hold on
                p(i,1:3) = plot_errorShade(x, targetERP_EP.mean, targetERP_EP.std, alpha, [0 0.60 1]);
                p(i,4:6) = plot_errorShade(x, targetERP.mean, targetERP.std, alpha, [1 0.20 0.60]);
                
                if j == 1
                    titStr = sprintf('%s\n%s\n%s', ['TARGET (', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')'], ...
                                    ['freq: ', num2str(parameters.oddballTask.targetFreq), ' Hz'], ...
                                    ['n = ', num2str(targetERP.n), ' / ', num2str(parameters.oddballTask.repeatsOfCycle)]);
                    tit(i) = title(titStr);
                else
                    tit(i) = title(['(', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')']);
                end
                
                if j == noOfChs
                    xLab(i) = xlabel('Time [ms]');
                else
                    xLab(i) = xlabel('');
                end
                yLab(i) = ylabel('\muV');

            i = 2;
            ind = (j-1)*cols+i;
                s(i) = subplot(rows, cols, ind);
                hold on
                p(i,1:3) = plot_errorShade(x, distracterERP_EP.mean, distracterERP_EP.std, alpha, [0 0.60 1]);
                p(i,4:6) = plot_errorShade(x, distracterERP.mean, distracterERP.std, alpha, [1 0.20 0.60]);
                
                if j == 1
                    titStr = sprintf('%s\n%s\n%s', ['DISTRACTER (', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')'], ...
                                    ['freq: ', num2str(parameters.oddballTask.distractFreq), ' Hz'], ...
                                    [' n = ', num2str(distracterERP.n), ' / ', num2str(parameters.oddballTask.repeatsOfCycle)]);
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
                p(i,1:3) = plot_errorShade(x, standardERP_EP.mean, standardERP_EP.std, alpha, [0 0.60 1]);
                p(i,4:6) = plot_errorShade(x, standardERP.mean, standardERP.std, alpha, [1 0.20 0.60]);
               
                if j == 1
                    titStr = sprintf('%s\n%s\n%s', ['STANDARD (', parameters.BioSemi.chName{ch + parameters.BioSemi.chOffset}, ')'], ...
                                    ['freq: ', num2str(parameters.oddballTask.standardFreq), ' Hz'], ...
                                    [' n = ', num2str(standardERP.n), ' / ', num2str(parameters.oddballTask.nrOfStandardsPerCycle*parameters.oddballTask.repeatsOfCycle)]);
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
                
                if j ==2
                    leg = legend([p(i,3) p(i,6)], 'EP Denoised', 'Filtered', 2, 'Location', 'NorthWest');
                        legend('boxoff')
                elseif j == 3
                    leg = legend([p(i,3) p(i,6)], ['ERP EP, scale: ', num2str(parameters.ep_den.scales_postStim)],...
                                 ['Filt. (',num2str(parameters.filter.bandPass_ERP_loFreq), '-', num2str(parameters.filter.bandPass_ERP_hiFreq), ' Hz)'], 2);
                    legend('boxoff')   
                        
                else
                    leg = legend('');
                        legend('hide')
                end

            % Style

                set(s, 'YLim', [-30 30])
                set(s, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-2)
                set([xLab yLab], 'FontName', style.fontName, 'FontSize', style.fontSizeBase-2, 'FontWeight', 'bold')
                set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-1, 'FontWeight', 'bold')
                set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-2, 'FontWeight', 'bold')
