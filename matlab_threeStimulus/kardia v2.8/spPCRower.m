 %% ---- spPCRower
function power=spPCRower(F, PSD, freq, parameters, callType)

    %% PETTERI: for LOMB-SCARGLE ESTIMATES with OUGP
    if strcmp(callType, 'OUGP')
        
        freqBinNames = fieldnames(parameters.heart.freqBins);
        indexBin = find(strcmp(freqBinNames, freq) == 1);
        noOfFreqBins = length(freqBinNames);
        
        % find indices corresponding the frequency indices
        for bin = indexBin : indexBin % noOfFreqBins       
            indexTemp = find(F == parameters.heart.freqBins.(freqBinNames{bin})(1));           
            if isempty(~indexTemp)
                %warning(['Start of the frequency bin ', num2str(parameters.heart.freqBins.(freqBinNames{bin})(1)), ' Hz not found for component "', freqBinNames{bin}, '"'])             
                %disp([freqBinNames{bin}, '  Check your frequency resolution settings (now ', num2str(parameters.heart.freqResolution), ' Hz or bin frequency ranges "parameters.heart.freqBins"'])
                [val1, fixed_ind1] = min(abs(F-parameters.heart.freqBins.(freqBinNames{bin})(1)));
                disp(['               .. ', freqBinNames{bin}, ' The closest start freq match is = ', num2str(val1), ' Hz (index = ', num2str(fixed_ind1), ')'])
                indexTemp = fixed_ind1;
            end
            freqIndices(bin,1) = indexTemp; % only one value should be found
            
            indexTemp = find(F == parameters.heart.freqBins.(freqBinNames{bin})(2));          
            if isempty(~indexTemp)
                %warning(['End of the frequency bin ', num2str(parameters.heart.freqBins.(freqBinNames{bin})(2)), ' Hz not found for component "', freqBinNames{bin}, '"'])
                %disp([freqBinNames{bin}, '  Check your frequency resolution settings (now ', num2str(parameters.heart.freqResolution), ' Hz or bin frequency ranges, "parameters.heart.freqBins"'])
                [val2, fixed_ind2] = min(abs(F - parameters.heart.freqBins.(freqBinNames{bin})(2)));
                disp(['               .. ', freqBinNames{bin}, ' The closest end freq match is = ', num2str(val2), ' Hz (index = ', num2str(fixed_ind2), ')'])
                indexTemp = fixed_ind2;
            end
            freqIndices(bin,2) = indexTemp; % only one value should be found
        end
        
            % The for maybe a bit weird here just for one component and was
            % mainly for all the components, but it was decided to match
            % the function with the output of Kardia so hopefully you are
            % okay with this
        
        % Integrate the PSD over the specified frequency range        
        power = parameters.heart.freqResolution * trapz(PSD(freqIndices(bin,1):freqIndices(bin,2)));
        

    %% KARDIA IMPLEMENTATION
    elseif strcmp(callType, 'kardia')

        % define frequency bands
        ulf = mean([0.0005 0.003]); % ultra-low frequency, http://dx.doi.org/10.1007/s11517-012-0928-2
        ulfStar = mean([0.002 0.01]); % ultra-low frequency STAR, http://dx.doi.org/10.1007/s11517-012-0928-2
        vlf=0.04; % very low frequency band
        lf=0.15; % low frequency band
        hf=0.4; % high frequency band

        % calculate number of points in the spectrum
        N=length(PSD);
        %calculate maximum frequency
        maxF=F(2)*N;

        if hf>F(end),
            hf=F(end);
            if lf>hf,
                lf=F(end-1);
                if vlf>lf,
                    vlf=F(end-2);
                    if ulfStar>vlf,
                        ulfStar=F(end-3);
                        if ulf>ulfStar,
                            ulf=F(end-4);
                        end
                    end
                end
            end
        end

        %calculate limiting points in each band
        index_ulf=round(ulf*N/maxF)+1;
        index_ulfStar=round(ulfStar*N/maxF)+1;
        index_vlf=round(vlf*N/maxF)+1;
        index_lf=round(lf*N/maxF)+1;
        index_hf=round(hf*N/maxF)+1;
        if index_hf>N,index_hf=N;end

        switch freq
            case {'total'}
                % calculate total energy (from 0 to hf) in ms^2
                total=F(2)*sum(PSD(1:index_hf-1));
                power=total;
            case {'ulf'}
                %calculate energy of very low frequencies (from 0 to uf2)
                ulf=F(2)*sum(PSD(1:index_ulf-1));
                power=ulf;
            case {'ulfStar'}
                %calculate energy of very low frequencies (from 0 to uf2)
                ulfStar=F(2)*sum(PSD(1:index_ulfStar-1));
                power=ulfStar;
            case {'vlf'}
                %calculate energy of very low frequencies (from 0 to vlf2)
                vlf=F(2)*sum(PSD(1:index_vlf-1));
                power=vlf;
            case {'lf'}
                %calculate energy of low frequencies (from vlf2 to lf2)
                lf=F(2)*sum(PSD(index_vlf:index_lf-1));
                power=lf;
            case {'hf'}
                %calculate energy of high frequencies (from lf2 to hf2)
                hf=F(2)*sum(PSD(index_lf:index_hf-1));
                power=hf;
            case {'nlf'}
                %calculate normalized low frequency
                lf=F(2)*sum(PSD(index_vlf:index_lf-1));
                hf=F(2)*sum(PSD(index_lf:index_hf-1));
                nlf=lf/(lf+hf);
                power=nlf;
            case {'nhf'}
                %calculate normalized low frequency
                lf=F(2)*sum(PSD(index_vlf:index_lf-1));
                hf=F(2)*sum(PSD(index_lf:index_hf-1));
                nhf=hf/(lf+hf);
                power=nhf;
            otherwise
                disp('Uknown frequency range selection')
                power=nan;
        end
    else
        callType
        error('Wrong callType, you have typo here?')
        
    end
end
    