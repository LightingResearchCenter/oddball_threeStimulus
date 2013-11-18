function epochOut = pre_removeBaseline_epochByEpoch(epochIn, j, parameters, handles)
        
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempBaselineRemoval.mat';
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

    debugEpochBaselineCorrection = 0;  

    % From ERPLAB:
    % http://erpinfo.org/erplab/erplab-documentation/documentation-archive-for-previous-versions/v1.x-documentation/erplab-manual/Epoching_Bins.htm
    % However, it does no harm to perform baseline correction at the epoching stage, 
    % and it may eliminate confusion if the epoched data are later exported to another system.  
    % Thus, we recommend that you perform baseline correction during the 
    % epoching process unless you have a good reason not to. 
    [rows,cols] = size(epochIn);
    x_in = (linspace(1,rows,rows))';  

    % baseline period indices, detrend based on the pre-stimulus EEG only:
    basIndex1 = parameters.oddballTask.baselineRemove_index1; 
    basIndex2 = parameters.oddballTask.baselineRemove_index2 * parameters.EEG.srate;
    baseline_x = (linspace(1,basIndex2,basIndex2))';
    
    epochOut_rmbas = zeros(rows,cols);
    epochOut_mean = zeros(rows,cols);

    % Remove Bassline
    for i = 1 : cols        
        baseline = (basIndex1:1:basIndex2);
        frames = length(epochIn(:,i));        
        [epochOut_rmbas(:,i), epochOut_mean(:,i)] = rmbase(epochIn(:,i)', frames, baseline);        
    end   

    epochOut = epochOut_rmbas;

    %% DEBUG (Alternative methods
    if debugEpochBaselineCorrection == 1

        epochOut_detr = zeros(rows,cols);
        epochOut_meanRemoved = zeros(rows,cols);

        % REMOVE BASELINE, CHANNEL BY CHANNEL
        for i = 1 : cols

            % JUST REMOVE THE MEAN OF PRE-STIMULUS (or any other time period
            % specified above)
            meanOfBaselinePeriod = nanmean(epochIn(basIndex1:basIndex2,i));

            % if epoch contains NaNs in baseline, then the whole epoch is an
            % artifact
            if isnan(meanOfBaselinePeriod)
                epochIn(:,i) = NaN;
            else
                epochOut_meanRemoved(:,i) = epochIn(:,i) - meanOfBaselinePeriod;
            end

            %{
            plot(x_in, epochIn(:,i), 'r', x_in, epochOut_meanRemoved(:,i), 'k')
            pause(2.0)
            %}

            % DETREND        

            % get the trend
            trendRemoved = detrend(epochIn(basIndex1:basIndex2,i), 'linear');
                % you could use 1st order polyfit as well, could be waster with
                % an analytic expression and no need to interpolate

            % subtract the removeTrend from input to get the actual trend as
            % the detrend cannot return it directly
            trend = epochIn(basIndex1:basIndex2,i) - trendRemoved;

            % now the trend has only a fraction (half to be exact if you use
            % the baseline before stimulus) of the original epoch datapoints so
            % we interpolate to the final length
            linearTrend = interp1(baseline_x, trend, x_in, 'linear', 'extrap');

            epochOut_detr(:,i) = epochIn(:,i) - linearTrend;

        end

        % RETURN
        % epochOut = epochOut_meanRemoved;
        % epochOut = epochOut_detr;

        % FROM ERPLAB:
        % ------------
        % http://erpinfo.org/erplab/erplab-documentation/documentation-archive-for-previous-versions/v1.x-documentation/erplab-manual/Epoching_Bins.htm

            % The Baseline correction option allows you to enable or disable baseline correction.  
            % If enabled, you can select the period that will be used for baseline correction. 
            % If you select Pre, the prestimulus baseline period will be used (this is the default).  
            % You could instead select Post to use the poststimulus period or Whole to select the entire epoch.  
            % Finally, you could select Custom and then provide two numbers that specify the beginning 
            % and end of the baseline period (e.g., -50 50 to use the interval from -50 ms to +50 ms).  
            % The baseline period must be entirely within the period of the epoch.  
            % For whatever period you select, the mean voltage over this period will be subtracted 
            % from the waveform for a each epoch (separately for each channel).

        epochOut_bp = zeros(rows,cols);    

            % Use short variable names for detrending filtering
            loFreq = 0.01; % parameters.epochERP.detrendingLowCut;
            hiFreq = 300; % parameters.epochERP.detrendingHighCut;
            N = 4; % parameters.epochERP.detrendingOrder;

            % Bandpass filter
            for i = 1 : cols
                epochOut_bp(:,i) = pre_bandbassFilter(epochIn(:,i), parameters.EEG.srate, [hiFreq, loFreq], N, N*10, handles);   
            end

        % Detrending
        for i = 1 : cols
            epochOut_detr(:,i) = detrend(epochIn(:,i), 'linear');   
        end

        %% Debug PLOT    
        debug_plotBaselineRemoval(epochOut_meanRemoved, epochOut_bp, epochOut_detr, epochOut_rmbas, epochOut_mean, epochIn, j, parameters, handles)            
                
    end

    function debug_plotBaselineRemoval(meanRemoved, bp, detr, rmbas, mean, In, j, parameters, handles)

        %{
        scrsz = get(0,'ScreenSize'); % get screen size for plotting

            fig = figure('Color', 'w');
                set(fig, 'Position', [0.05*scrsz(3) 0.20*scrsz(4) 0.92*scrsz(3) 0.62*scrsz(4)])
        %}

        t = ((1 : 1 : length(meanRemoved))' / parameters.EEG.srate) - parameters.oddballTask.ERP_baseline;

        subplot(2,1,1)
            hold on
            plot(t, meanRemoved(:,1:1), 'r')
            plot(t, bp(:,1:1), 'g')
            plot(t, detr(:,1:1), 'b')
            plot(t, rmbas(:,1:1), 'k')                
                l(1) = legend('meanRemoved', 'BandPass', 'Detrend', 'rmbase');
                    legend('boxoff')
                    xlabel('Time [s]'); ylabel('Amplitude \muV');
                title(['Epoch #', num2str(j)])
                line([min(t) max(t)], [0 0], 'Color', 'k')
            hold off

        subplot(2,1,2)
            hold on
            plot(t, In(:,1:1), 'r')
            plot(t, mean(:,1:1), 'm') 
            hold off
                l(2) = legend('In', 'mean');
                    legend('boxoff')
                    xlabel('Time [s]'); ylabel('Amplitude \muV');

        set(l, 'Location', 'NorthEastOutside')

        if handles.figureOut.debugON == 1    
            i = 1;
            try
                i = 2;
                drawnow
                dateStr = plot_getDateString(); % get current date as string                               
                %fileNameOut = ['debug_baselineRemoval_', strrep(handles.inputFile, '.pdf', ''), '_', dateStr];
                fileNameOut = ['debug_baselineRemoval_', dateStr];
                % disp([' ... saving figure to disk (', fileNameOut, '.png]'])
                % export_fig(fullfile(handles.path.debugPreprocessing, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)                                    
            catch err
                err
            end
        end
