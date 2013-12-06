 function [f, amplit, PSD, amplit_matrix, PSD_matrix] = analyze_powerSpectrum(EEGchannel, Fs, windowType, r, segmentLength, nfft, freqRange, nOverlap, fromWhereCalled)
       
    overlapInSamples = segmentLength * (1 - (nOverlap/100));    
    nrOfWindows = 1 + ( (length(EEGchannel)-segmentLength) / overlapInSamples);  % e.g. 1 + (333824 - (2*2048)) / (1*2048) = 1 + 148 = 149

    x = 1 : 1 : segmentLength/2 + 1;
    f_in = freqRange;
    f = f_in;
    
    if strcmp(windowType, 'Tukey')
        window = (tukeywin(segmentLength,r)); % Tukey window with r=0.10 is 10% cosine window      
    else
        windowType
        error('Only Tukey implemented at the moment')
    end

    % Define manually the segments for the pwelch routine, as our EEG
    % data most likely contain NaNs for artifacts, the "automatic
    % overlapping" of pwelch will give us NaN values

    % preallocate
    PSD_matrix = zeros(length(x), nrOfWindows);
    amplit_matrix = zeros(length(x), nrOfWindows);       

    if strcmp(fromWhereCalled, 'eegTimeSeries') == 1
        fprintf('        .. ')
    end
    
    for i = 1 : nrOfWindows

        if rem(i,floor(nrOfWindows/6)) == 0 || i == nrOfWindows
            if strcmp(fromWhereCalled, 'eegTimeSeries') == 1
                fprintf('%s%s', num2str(100*(i/nrOfWindows),3), '% ')
                if i == nrOfWindows
                    fprintf('%s\n', ' ')
                end
            end
        end

        % get indices for this Window
        i1 = (i-1) * overlapInSamples + 1;
        i2 = i1 + segmentLength - 1;
        
        % will return empty power spectrum if you have some NaN values in
        % the input
        % numberOfNaNs = length(EEGchannel(isnan(EEGchannel(i1:i2))));
        
        % Mean-square power spectrum                        
        [PSD_matrix(:,i), f1] = pwelch(EEGchannel(i1:i2),segmentLength,nOverlap,nfft,Fs, 'PSD');     
        
        
        % Power spectral density
        [amplit_matrixTemp, f2] = pwelch(EEGchannel(i1:i2),segmentLength,nOverlap,nfft,Fs,'centered','power');
        
        % DEBUG
        %{
        subplot(2,1,1)
        plot(EEGchannel(i1:i2)); title(num2str(i)); drawnow
        subplot(2,1,2)
        plot(amplit_matrixTemp); title(num2str(i)); drawnow
        pause(0.5)
        %}

            % now the amplit_matrixTemp is two-sided amplitude estimate,
            % trim half a way
            iStart = segmentLength /2;
            amplit_matrix(:,i) = amplit_matrixTemp(iStart:end);

    end  

    
    
    PSD.mean = nanmean(PSD_matrix,2);
    PSD.SD = nanstd(PSD_matrix,1,2);
    amplit.mean =  nanmean(amplit_matrix,2);
    amplit.SD =  nanstd(amplit_matrix,1,2);
    
    try       
        f = f1;
    catch err
        err
        whos
        error('No f1 from PSD computation? The epoch windowing problem, fix if needed')
    end
        
       

