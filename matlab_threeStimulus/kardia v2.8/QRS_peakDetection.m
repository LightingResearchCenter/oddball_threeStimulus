function [rPeakTimes, rPeakAmplitudes] = QRS_peakDetection(ecg_data, fs)

    % from: http://matlabz.blogspot.com/2011/04/contents-cancellation-dc-drift-and.html
    
    % QRS Detection Example
    % shows the effect of each filter according to Pan-Tompkins algorithm.
    % Note that, the decision  algorithm is different then the mentioned algorithm.
    % by Faruk UYSAL


    %ecg_data = load('ecg3.dat'); % load the ECG signal from the file
    %fs = 200;              % Sampling rate
    N = length(ecg_data);       % Signal length
    t = [0:N-1]/fs;        % time index
    debugPlot = 1;


    %% DEBUG PLOT
    if debugPlot == 1
        fig = figure('name', 'QRS: Pan-Tompkins', 'Color', 'w');
        scrsz = get(0,'ScreenSize'); % get screen size for plotting
            set(fig, 'Position', [0.3*scrsz(3) 0.05*scrsz(4) 0.55*scrsz(3) 0.90*scrsz(4)])
            rows = 5;
            cols = 2;
        
        spIndex = 1;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot(t,ecg_data)
        xlabel('second');ylabel('Volts');title('Input ECG Signal')
    end
        

    %% CANCELLATION DC DRIFT AND NORMALIZATION

    ecg_data = ecg_data - mean (ecg_data );    % cancel DC conponents    
    ecg_data = detrend(ecg_data); % PETTERI, added detrending
    ecg_data = ecg_data/ max( abs(ecg_data )); % normalize to one

    if debugPlot == 1
        spIndex = 2;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot(t,ecg_data)
        xlabel('second');ylabel('Volts');title(' ECG Signal after cancellation DC drift and normalization')
    end
    

    %% LOW PASS FILTERING

    % LPF (1-z^-6)^2/(1-z^-1)^2
    b=[1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a=[1 -2 1];


    h_LP=filter(b,a,[1 zeros(1,12)]); % transfer function of LPF

    x2 = conv (ecg_data ,h_LP);
    %x2 = x2 (6+[1: N]); %cancle delay
    x2 = x2/ max( abs(x2 )); % normalize , for convenience .

    if debugPlot == 1
        spIndex = 3;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x2)-1]/fs,x2)
        xlabel('second');ylabel('Volts');title(' ECG Signal after LPF')
        xlim([0 max(t)])
    end

    %% HIGH PASS FILTERING

    % HPF = Allpass-(Lowpass) = z^-16-[(1-z^-32)/(1-z^-1)]
    b = [-1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 32 -32 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1];
    a = [1 -1];

    h_HP=filter(b,a,[1 zeros(1,32)]); % impulse response iof HPF

    x3 = conv (x2 ,h_HP);
    %x3 = x3 (16+[1: N]); %cancle delay
    x3 = x3/ max( abs(x3 ));

    if debugPlot == 1
        spIndex = 4;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x3)-1]/fs,x3)
        xlabel('second');ylabel('Volts');title(' ECG Signal after HPF')
        xlim([0 max(t)])
    end
    

    %% DERIVATIVE FILTER

    % Make impulse response
    h = [-1 -2 0 2 1]/8;
    % Apply filter
    x4 = conv (x3 ,h);
    x4 = x4 (2+[1: N]);
    x4 = x4/ max( abs(x4 ));

    % trim the start off (or turn to NaNs)
    x4(1:1000) = NaN;
    %x4 = x4(60:end);
    N = length(x4);
    
    if debugPlot == 1
        spIndex = 5;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x4)-1]/fs,x4)
        xlabel('second');ylabel('Volts');title(' ECG Signal after Derivative')
    end
    
    
    

    %% SQUARING
    x5 = x4 .^2;
    x5 = x5/ max( abs(x5 ));
    
    if debugPlot == 1
        spIndex = 6;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x5)-1]/fs,x5)
        xlabel('second');ylabel('Volts');title(' ECG Signal Squareing')
    end

    %% MOVING WINDOW INTEGRATION

    % Make impulse response
    h = ones (1 ,31)/31;
    Delay = 15; % Delay in samples

    % Apply filter
    x6 = conv (x5 ,h);
    x6 = x6 (15+[1: N]);
    x6 = x6/ max( abs(x6 ));

    if debugPlot == 1
        spIndex = 7;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot([0:length(x6)-1]/fs,x6)
        xlabel('second');ylabel('Volts');title(' ECG Signal after Averaging')
    end
    
    %% FIND QRS POINTS WHICH IT IS DIFFERENT THAN PAN-TOMPKINS ALGORITHM
    max_h = max(x6);
    thresh = mean (x6 );
    poss_reg =(x6>thresh*max_h)';

    if debugPlot == 1
        spIndex = 8;
        sp(spIndex) = subplot(rows,cols,spIndex);
        hold on
        plot (t, ecg_data/max(ecg_data))
        box on
        xlabel('second');ylabel('Integrated')
        xlim([1 3])
    end

    left = find(diff([0 poss_reg])==1);
    right = find(diff([poss_reg 0])==-1);    
   
    left=left-(6+16); % cancle delay because of LP and HP
    right=right-(6+16); % cancle delay because of LP and HP

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

    for i=1:length(left)

        try
            [R_value(i) R_loc(i)] = max( ecg_data(left(i):right(i)) );            
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
        
        spIndex = 9;
        sp(spIndex) = subplot(rows,cols,spIndex);
        title('ECG Signal with R points');
        plot (t,ecg_data/max(ecg_data) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o');
        legend('ECG','R','S','Q');
        legend('boxoff')
        
        spIndex = 10;
        sp(spIndex) = subplot(rows,cols,spIndex);
        plot (t,ecg_data/max(ecg_data) , t(R_loc) ,R_value , 'r^', t(S_loc) ,S_value, '*',t(Q_loc) , Q_value, 'o');
        xlim([1 3])
    end
    
    %% RETURN R peaks
    rPeakTimes = t(R_loc);
    rPeakAmplitudes = R_value;
    
        % not really a good implementation with high sample rates (e.g. 4,096 Hz).. returns the actual R peaks but
        % they come with additional peaks :S (Petteri Teikari)
        
        
