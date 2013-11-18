function [RMS_alphaFixed, RMS_alphaIAF] = pre_baselinePowerAnalysis(epoch, baselineCorr, alpha, parameters, handles)
            
    % Barry RJ, Kirkaikul S, Hodder D. 2000. 
    % EEG alpha activity and the ERP to target stimuli in an auditory oddball paradigm. 
    % International Journal of Psychophysiology 39:39–50. 
    % http://dx.doi.org/10.1016/S0167-8760(00)00114-8.
    segmentLength = parameters.oddballTask.preERP_power_segmentLength*parameters.EEG.srate;                    
    nfft = segmentLength;
    freqRange = 0 : 0.1 : parameters.filter.bandPass_hiFreq;
    nOverlap = parameters.oddballTask.preERP_power_nOverlap;                    
    preStimulusData = epoch(1:baselineCorr);

    % use the same subfunction as with the general
    % timeSeriesEEG analysis
    [f, amplit, PSD, ~, ~] = analyze_powerSpectrum(preStimulusData, parameters.EEG.srate, 'Tukey', ...
                                                parameters.oddballTask.preERP_power_tukeyWindowR, segmentLength, ...
                                                nfft, freqRange, nOverlap, 'preStimulus');

    % Construct the 'alpha power spectrum' as now the RMS                                         

    % get scalar RMS amplitude from the power spectra
    % calculated, now the frequency resolution is not great
    % as the segment is only 0.5 seconds 
    [RMS_alphaFixed, RMS_alphaIAF] = analyze_getBandAmplitude(f, amplit.mean, amplit.SD, parameters.powerAnalysis.alphaRange, 'Power', ...
                                     alpha.amplit_gravity, parameters.oddballTask.preERP_IAF_range, handles);

                                 % gravity recommended over
                                 % peak by Klimesch 1999,
                                 % you could add a switch
                                 % here if you wanna play
                                 % with peak