function [list_properties, activData, blinkData] = component_properties_mod(EEG, blink_chans, lpf_band)
    
    %{
    dir = cd;
    if nargin == 0
        load(fullfile(dir, 'debug.mat'))
    else
        save(fullfile(dir, 'debug.mat'))
    end    
    %}


% Copyright (C) 2010 Hugh Nolan, Robert Whelan and Richard Reilly, Trinity College Dublin,
% Ireland
% nolanhu@tcd.ie, robert.whelan@tcd.ie
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

dataIn = EEG.data;
activ= EEG.icaact;

list_properties = [];
%
if isempty(EEG.icaweights)
    fprintf('No ICA data.\n');
    return;
end

if ~exist('lpf_band','var') || length(lpf_band)~=2 || ~any(lpf_band)
    ignore_lpf=1;
else
    ignore_lpf=0;
end

delete_activations_after=0;
if ~isfield(EEG,'icaact') || isempty(EEG.icaact)
    delete_activations_after=1;
    EEG.icaact = eeg_getica(EEG);
end

for u = 1:size(EEG.icaact,1)
    [spectra(u,:) freqs] = pwelch(EEG.icaact(u,:),[],[],(EEG.srate),EEG.srate);
end

list_properties = zeros(size(EEG.icaact,1),5); % This 5 corresponds to number of measurements made.

[noOfICA_PCA_channels, noOfSamples] = size(EEG.icaact)
for u=1:noOfICA_PCA_channels

    measure = 1;
    %% TEMPORAL PROPERTIES

        % 1 Median gradient value, for high frequency stuff
        list_properties(u,measure) = median(diff(EEG.icaact(u,:)));
        measure = measure + 1;

        % 2 Mean slope around the LPF band (spectral)
        if ignore_lpf
            list_properties(u,measure) = 0;
        else
            try
                meanSlope = mean(diff(10*log10(spectra(u,find(freqs>=lpf_band(1),1):find(freqs<=lpf_band(2),1,'last')))));
                list_properties(u,measure) = meanSlope;
            catch err
                %err
                if isempty(meanSlope)
                    disp('      ... freqVector empty, no slope done, adjust low-pass band')
                    ignore_lpf = 1;
                end
            end
        end
        measure = measure + 1;

    %% SPATIAL PROPERTIES

        % 3 Kurtosis of spatial map (if v peaky, i.e. one or two points high
        % and everywhere else low, then it's probably noise on a single
        % channel)
        list_properties(u,measure) = kurt(EEG.icawinv(:,u));
        measure = measure + 1;

    %% OTHER PROPERTIES

        % 4 Hurst exponent
        list_properties(u,measure) = hurst_exponent(EEG.icaact(u,:));
        measure = measure + 1;

    %% 10 Eyeblink correlations
    memWarningFlag = 0;
    
        if (exist('blink_chans','var') && ~isempty(blink_chans))
            for v = 1:length(blink_chans)
                if ~(max(EEG.data(blink_chans(v),:))==0 && min(EEG.data(blink_chans(v),:))==0);                
                    blinkData = EEG.data(blink_chans(v),:);
                    activData = EEG.icaact(u,:);                
                    try
                        % plot(activData, blinkData, 'o')
                        f = corrcoef(EEG.icaact(u,:), EEG.data(blink_chans(v),:));
                        
                    catch err

                        if strcmp(err.identifier, 'MATLAB:nomem')

                           blinkLength = length(EEG.data(blink_chans(v),:));
                           activLength = length(EEG.icaact(u,:));

                           if memWarningFlag == 0

                               mem_Mb_needed = (activLength * activLength * 4) / 1024 / 1024;
                               warnStr = ['   .. Not enough MEMORY available: ', num2str(blinkLength), ' x ', num2str(activLength), ...
                                        ' x 4 bytes (single) = ', num2str(mem_Mb_needed), ' MB needed'];
                               disp(warnStr)
                               memWarningFlag = 1;

                               % for Windows, you could use the memory
                               % memory
                               % http://www.mathworks.com/help/matlab/matlab_prog/resolving-out-of-memory-errors.html

                           end

                           memTarget = 2048; % MB
                           downsamplingNeeded = mem_Mb_needed / memTarget;
                           downsamplingNeeded = 2^nextpow2(downsamplingNeeded); % find the next power of 2
                           x = linspace(1, blinkLength, blinkLength); 
                           x_i = linspace(1, blinkLength, blinkLength/downsamplingNeeded);

                           disp(['      ... downsampling (1/', num2str(downsamplingNeeded), ') the blink channel and the ICA activation'])
                           blinkData = interp1(x, EEG.data(blink_chans(v),:), x_i);
                           activData = interp1(x, EEG.icaact(u,:), x_i);       

                           % try again to compute the correlation coefficient
                           f = corrcoef(activData, blinkData);

                        else % unidentified error
                            %'MATLAB:samelen'                       
                            whos
                        end
                    end

                    if sum(sum(isnan(f))) == length(f)*length(f)
                        warning('You get werid behavior here if you have multiple corrcoef.m in your Matlab path!')
                        which -all corrcoef % https://groups.google.com/forum/#!topic/comp.soft-sys.matlab/cC77XLgs7_4
                        disp('    you will get a conflicting function annoyingly from BioSig package, so you have to remove the /NaN/inst/ from path')                
                    end                
                    x(v) = abs(f(1,2));
                else
                    x(v) = v;
                end
            end
            list_properties(u,measure) = max(x);
            measure = measure + 1;
        end
end

for u = 1:size(list_properties,2)
    list_properties(isnan(list_properties(:,u)),u)=nanmean(list_properties(:,u));
    list_properties(:,u) = list_properties(:,u) - median(list_properties(:,u));
end

if delete_activations_after
    EEG.icaact=[];
end