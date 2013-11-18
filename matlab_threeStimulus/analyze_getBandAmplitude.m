function [fixedBand, indivBand] = analyze_getBandAmplitude(f, ps_mean, ps_SD, freqRange, dataType, IAF, IAF_range, handles)

    debugMatFileName = 'tempBandAmplitude.mat';
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

    % condition the IAF, so that it is rounded to one decimal precision
    precision = 10;
    IAF_range(1) = round(precision*(IAF_range(1) + IAF)) / precision;
    IAF_range(2) = round(precision*(IAF_range(2) + IAF)) / precision;
    
    % Get indices, note that now it is possible that exact frequency
    % matches are not found, so we want the closest index then
    [val1,fixed_ind1] = min(abs(f-freqRange(1)));
    [val2,fixed_ind2] = min(abs(f-freqRange(2)));
    [val3,indiv_ind1] = min(abs(f-IAF_range(1)));
    [val4,indiv_ind2] = min(abs(f-IAF_range(2)));
    
    % Get the powers / PSD of the specified freq band
        
        % for fixed frequency range
        amplit_fixedBand = ps_mean(fixed_ind1:fixed_ind2);
        f_fixed = f(fixed_ind1:fixed_ind2);
    
        % for the individual frequency range
        amplit_indivBand = ps_mean(indiv_ind1:indiv_ind2);
        f_indiv = f(indiv_ind1:indiv_ind2);
        
        freqRes = f_fixed(2) - f_fixed(1);    
    
    if strcmp(dataType, 'Power')
        
        % Calculate the Amplititude
        fixedBand = sqrt(sum(amplit_fixedBand .^ 2));
        indivBand = sqrt(sum(amplit_indivBand .^ 2));
        
    elseif strcmp(dataType, 'PSD') 
        
        % Calculate the Power, maybe change the variable names later, but
        % outside this function the input data is going to be PSD and not
        % amplitude
        fixedBand = freqRes * trapz(amplit_fixedBand); 
        indivBand = freqRes * trapz(amplit_indivBand);
        
    else err
        err
        error(['Error with your dataType'])
    end
        
    % Debug
        if handles.flags.showDebugMessages == 1
            disp(['        IAF Frequency: ', num2str(IAF,4), ' Hz'])
            disp(['          Fixed range (match): ', num2str(f(fixed_ind1)), '-', num2str(f(fixed_ind2)), ' Hz'])
            disp(['            Fixed RMS: ', num2str(fixedBand), ' uV'])
            disp(['              IAF range (match): ', num2str(f(indiv_ind1)), '-', num2str(f(indiv_ind2)), ' Hz'])
            disp(['                IAF RMS: ', num2str(indivBand), ' uV'])
        end

        if handles.flags.showDebugPlots == 1
            plot(f_fixed, amplit_fixedBand, 'b', f_indiv, amplit_indivBand+1, 'r')
                legend('Fixed', 'IAF', 2, 'Location', 'Best'); legend('boxoff')
                xlabel('Hz'); ylabel('\muV');
                titStr = sprintf('%s\n%s', ['input freq resolution: ', num2str(f(2)-f(1)), ' Hz'], 'vertical displacement for visualization');
                title(titStr)
        end
