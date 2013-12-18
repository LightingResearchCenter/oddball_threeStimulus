function [rPeakTimes, rPeakAmplitudes, forDebugPlot_QRS] = QRS_peakDetection(ecg_data, Fs, handles)

    % from: http://matlabz.blogspot.com/2011/04/contents-cancellation-dc-drift-and.html
    
    % QRS Detection Example
    % shows the effect of each filter according to Pan-Tompkins algorithm.
    % Note that, the decision  algorithm is different then the mentioned algorithm.
    % by Faruk UYSAL

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempQRS_Faruk.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
            % for debugging, use shorter vector
            %ecg_data = ecg_data(1:4096*60);
        else
            if handles.flags.saveDebugMATs == 1
                path = handles.path;
                save('debugPath.mat', 'path')
                save(fullfile(path.debugMATs, debugMatFileName))            
            end
        end 
    end

    
    
    %ecg_data = load('ecg3.dat'); % load the ECG signal from the file
    %Fs = 200;              % Sampling rate
    N = length(ecg_data);       % Signal length
    t = [0:N-1]/Fs;        % time index
    debugPlot = 1;
    zoomWindow = [5 6.6];
    whos


    %% DEBUG PLOT
    if debugPlot == 1
        fig = figure('name', 'QRS: Pan-Tompkins', 'Color', 'w');
        scrsz = get(0,'ScreenSize'); % get screen size for plotting
            set(fig, 'Position', [0.03*scrsz(3) 0.05*scrsz(4) 0.55*scrsz(3) 0.90*scrsz(4)])
            rows = 4;
            cols = 3;
        
        spIndex = 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot(t,ecg_data,'Color',[0 0.6 1])
        xlabel('second');ylabel('\muVolts');title('Input ECG Signal')
        xlim(zoomWindow)
    end
        

    %% CANCELLATION DC DRIFT AND NORMALIZATION

    ecg_data = ecg_data - mean (ecg_data );    % cancel DC components    
    ecg_data = detrend(ecg_data); % PETTERI, added detrending
    ecg_data = ecg_data/ max( abs(ecg_data )); % normalize to one

    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot(t,ecg_data,'Color',[0 0.6 1])
        xlabel('second');ylabel('\muVolts');title(' ECG Signal after cancellation DC drift and normalization')
        xlim(zoomWindow)
    end
    

    %% LOW PASS FILTERING & HIGH PASS FILTERING

        % PETTERI: Combine the steps together from original implementation
        % that used LP + HP with some hard-coded delay correction, now we
        % use ZERO-PHASE filter to preserve the temporal characteristics of
        % the ECG waveform
    
        cutOffs = [15 5];
        filterOrder = 10; % conservative 
        try 
            x2 = pre_bandbassFilter(ecg_data, Fs, cutOffs, filterOrder, filterOrder, handles);
        catch err
            err
            filterOrder = filterOrder / 2;
            x2 = pre_bandbassFilter(ecg_data, Fs, cutOffs, filterOrder, filterOrder, handles);
            error('Bandpass filter parameters too "extreme"')
        end
        
        % NOTE! We already have once bandpass filtered the ECG data in
        % pre_denoiseECG.m but then the frequency range was 2.3 - 30 Hz,
        % and now we are using a lot tighter 5-11 Hz bandpass, see e.g.
        % http://ocw.utm.my/file.php/38/SEB4223/07_ECG_Analysis_1_-_QRS_Detection.ppt%20%5BCompatibility%20Mode%5D.pdf
        % Slides by Dr. Malarvili Balakrishnan on Pan-Tompkins
        
    %x2 = x2/ max( abs(x2 )); % normalize , for convenience .

    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x2)-1]/Fs,x2,'Color',[0 0.6 1])
        xlabel('second');ylabel('\muVolts');title([' ECG after Bandpass [', num2str(cutOffs(2)), '-', num2str(cutOffs(1)), ' Hz]'])
        xlim(zoomWindow)
    end


    %% DERIVATIVE FILTER

    x3 = x2; % now we modified to code so the x3 came from the separate HIGHPASS FILTER
    
    % Instead of using the original derivative filter in the implementation    
        h = [-1 -2 0 2 1]/8;
        
    % We could use the Smooth Differentiation by Jianwen Luo, not to
    % accentuate the noise too much
    % http://www.mathworks.com/matlabcentral/fileexchange/6170-smooth-differentiation 
        %diff_filterLength = 32*16;
        %h = smooth_diff(diff_filterLength);
        
    % Apply filter
    x4 = conv (x3 ,h);
    x4 = x4 (2+[1: N]);
    x4 = x4/ max( abs(x4 ));

    % trim the start off (or turn to NaNs)
    x4(1:1000) = NaN;
    %x4 = x4(60:end);
    N = length(x4);
    
    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x4)-1]/Fs,x4,'Color',[0 0.6 1])
        xlabel('second');ylabel('\muVolts');title(' ECG Signal after Derivative')
        xlim(zoomWindow)
    end
    drawnow
 
    

    %% SQUARING
    x5 = x4 .^2;
    x5 = x5/ max( abs(x5 ));
    
    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x5)-1]/Fs,x5,'Color',[0 0.6 1])
        xlabel('second');ylabel('\muVolts');title(' ECG Signal Squareing')
        xlim(zoomWindow)
    end

    %% MOVING WINDOW INTEGRATION   

    % Apply filter (original convolution)    
    filtLength = 31;
    h = ones (1,filtLength)/31; % Make impulse response
    Delay = floor(filtLength/2); % Delay in samples
    
    x6_init = conv (x5,h);    
    x6 = x6_init(Delay+[1: N]);         
    x6 = x6/ max( abs(x6));
   
    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        hold on
        multip = 20;
        plot(x5 / max(x5), 'r')
        plot(x6_init / max(x6_init), 'g')
        plot(x6 / max(x6), 'k')
        hold off
        %leg = legend(['Sq., n=', num2str(length(x5))], ['x6init, n=', num2str(length(x6))], ['x6, n=', num2str(length(x6))]);
        leg = legend(['Sq.'], ['x6init'], ['x6']);
            set(leg, 'FontSize', 6, 'Location', 'Best');
            legend('boxoff'); title('Average')            
            xlim([Delay*multip*3 Delay*multip*8])    
            
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x6)-1]/Fs,x6,'Color',[0 0.6 1])
        xlabel('second');ylabel('\muVolts');title(' ECG Signal after Averaging')
        xlim(zoomWindow)
    end
    
    %% FIND QRS POINTS WHICH IT IS DIFFERENT THAN PAN-TOMPKINS ALGORITHM
    max_h = nanmax(x6);
    thresh = nanmean(x6);
    poss_reg =(x6>thresh*max_h)';
    
    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        hold on
        plot (t, ecg_data/max(ecg_data),'Color',[0 0.6 1])
        box on
        xlabel('second'); ylabel('\muVolts'); title('ECG')
        xlim(zoomWindow)
    end

    if debugPlot == 1
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        hold on
        plot(t, x6 - (thresh*max_h),'Color',[0 0.6 1])
        plot (t, ecg_data/max(ecg_data),'Color',[0.8 0.2 1])
        title(['max_h = ', num2str(max_h), ', thresh = ', num2str(thresh)])
        box on
        xlim(zoomWindow)
            leg = legend('Thresh*max_h', 'ECG', 'Location', 'best'); 
                legend('boxoff'); set(leg, 'FontSize', 7)
        drawnow
    end
    
    left = find(diff([0 poss_reg])==1);
    right = find(diff([poss_reg 0])==-1);    
   
    %left=left-(6+16); % cancle delay because of LP and HP
    %right=right-(6+16); % cancle delay because of LP and HP

        % PETTERI: correct for negative values
        leftNegative = left < 0;
        numberOfLeftNegatives = sum(leftNegative);
        leftTemp = left(~leftNegative);    
        rightTemp = right(~leftNegative);

        rightNegative = rightTemp < 0;
        numberOfRightNegatives = sum(rightNegative);
        leftTemp = leftTemp(~rightNegative);
        rightTemp = rightTemp(~rightNegative);
        
        left = leftTemp;
        right = rightTemp;

    if debugPlot == 1
        
        spIndex = 10;
        sp(spIndex) = subplot(rows,cols,spIndex);
        rightVector = zeros(length(t),1); rightVector(right) = 1;
        leftVector = zeros(length(t),1); leftVector(left) = 1;
        hold on
        plot (t, ecg_data/max(ecg_data),'Color',[0.8 0.2 1])
        plot(t, leftVector, 'Color',[0 0.6 1])
        plot(t, rightVector, 'Color',[0.3 0.9 0.4])        
        title('Search Window')
        box on
        xlim(zoomWindow)
            leg = legend('ECG', 'Left', 'Right', 'Location', 'best'); 
                legend('boxoff'); set(leg, 'FontSize', 7)
        
    end
        
    for i=1:length(left)

        try
            [R_value(i) R_loc(i)] = max(ecg_data(left(i):right(i)) );            
        catch err
            err
            left(i)
            right(i)
        end
        
        R_loc(i) = R_loc(i)-1+left(i); % add offset

        [Q_value(i) Q_loc(i)] = min( ecg_data(left(i):R_loc(i)) );
        Q_loc(i) = Q_loc(i)-1+left(i); % add offset

        [S_value(i) S_loc(i)] = min( ecg_data(left(i):right(i)) );
        S_loc(i) = S_loc(i)-1+left(i); % add offset

    end

    % there is no selective wave
    try 
        Q_loc=Q_loc(Q_loc~=0);
    catch err
        % err
        disp('            no Q_loc')
    end
    
    try 
        R_loc=R_loc(R_loc~=0);
    catch err
        % err
        rPeakTimes = NaN;
        rPeakAmplitudes = NaN;        
        return
    end
    S_loc=S_loc(S_loc~=0);

    %% FINAL PLOT
    if debugPlot == 1
        
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        title('ECG Signal with R points');
        p1 = plot(t, ecg_data, t(R_loc) ,R_value, 'r^', t(S_loc) ,S_value, 'b*', t(Q_loc) , Q_value, 'go');
            set(p1(1), 'Color', [.7 .7 .7])        
        xlim([0 max(t)])
        
        spIndex = spIndex + 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        p2 = plot(t,ecg_data, t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, 'b*',t(Q_loc) , Q_value, 'go');
            set(p2(1), 'Color', [.7 .7 .7])
        xlim(zoomWindow)
        
        leg = legend('ECG','R','S','Q');
            set(leg,'Position',[0.902867965367965 0.183559782608696 0.0871212121212121 0.0944293478260869], 'FontSize', 7);
            legend('boxoff')
            drawnow
        
        % Auto-SAVE
        try
            if handles.figureOut.ON == 1                     
                dateStr = plot_getDateString(); % get current date as string          
                fileNameOut = sprintf('%s%s', 'debug_QRS_panTomp_', strrep(handles.inputFile, '.bdf', ''), '.png');
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
            end
        catch err
            err
            str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                          'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
            error(str)
        end        
        
    end
    
    %% DO FINAL OUTLIER EXCLUSION
    
        % We can use the same algorithm used after this to reject outlier RR
        % intervals, and it is actually more beneficial to do it here. In
        % practice the algorithm may find R peaks that have too small
        % amplitudes and we could check for deviating R peak amplitudes
        rPeakTimes = t(R_loc);  
        handles.parameters.heart.ectopicMedianParam = 12;
        handles.parameters.heart.ectopicOutlierMethod = 'median'; % 'percent', 'median', 'sd'    
        handles.parameters.heart.outlierReplaceMethod = 'remove'; % 'mean' / 'median' / 'cubic' / 'remove'
        
        disp(['           * Rejecting outliers from R peaks (only R, not Q or S)'])
        [rPeakAmplitudes, rPeakTimes] = pre_correctHeartRatePeriodForOutliers([], rPeakTimes, R_value, 'rPeak', handles);                   
        noOfPeaksRejected = length(t(R_loc)) - length(rPeakTimes);
        
        if debugPlot == 1
            
            axes(sp(spIndex-1))
                hold on
                plot(rPeakTimes, rPeakAmplitudes, 'k^')
                hold off
                title(['QRS, ', num2str(noOfPeaksRejected), ' x R rejected'])
            
            axes(sp(spIndex))
                hold on
                plot(rPeakTimes, rPeakAmplitudes, 'k^')
                hold off
                title('QRS Zoom')
            
                leg = legend([{'ECG'}, {'R w outl'}, {'S'}, {'Q'}, {'R'}]);
                    set(leg,'Position',[0.902867965367965 0.183559782608696 0.0871212121212121 0.0944293478260869], 'FontSize', 7);
                    legend('boxoff')
                    drawnow
            
            
        end
        % not really a good implementation with high sample rates (e.g. 4,096 Hz).. returns the actual R peaks but
        % they come with additional peaks :S (Petteri Teikari)
           
    %% RETURN DEBUG for plotting
    
        % So that you can plot all the heart parameters to a single plot
        % easier visual inspection afterwards 
        forDebugPlot_QRS.t = t;
        forDebugPlot_QRS.ECG = ecg_data;
        forDebugPlot_QRS.t_Rraw = t(R_loc);
        forDebugPlot_QRS.Rraw = R_value;
        forDebugPlot_QRS.t_Q = t(Q_loc);
        forDebugPlot_QRS.Q = Q_value;
        forDebugPlot_QRS.t_S = t(S_loc);
        forDebugPlot_QRS.S = S_value;
        forDebugPlot_QRS.t_R = rPeakTimes;
        forDebugPlot_QRS.R = rPeakAmplitudes;
        
        
