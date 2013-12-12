function [WTout, nonNormFreqIndex] = analyze_normalizeWaveletSpectrum(WT, WT_w_COI, t, f, baselineLimits, timeResolutionDivider, parameters, handles)
    
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempTFnorm.mat';
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

    timeRes = t(2) - t(1);
    freqRes = f(2) - f(1);   
    sampleRate = parameters.EEG.srate / timeResolutionDivider;
    %baselineLimits
    debugOn = 1;
    
    
    %% COMPUTATION   
        
        % find indices corresponding to the baseline time
        indicesTrim = analyze_getTFtrimIndices(t, baselineLimits, [], [], parameters);

        % get the baseline mean (one mean for each frequency)    
        meansPerFreq = nanmean(WT(:,indicesTrim(1):indicesTrim(2)),2); 
        baselineVector = meansPerFreq;

            nonNormFreqIndex = find(isnan(meansPerFreq) == 1, 1, 'last');        

        % a bit akward as now these values are not normalized, low frequencies
        % in other words if you don't normalize all the way to the stimulus
        % onset (e.g. between -500 and -100 ms)
        %baselineVector(isnan(baselineVector)) = 1;     

        % replicate the vector to a matrix so that the matrix dimensions match
        baselineVector = repmat(baselineVector, 1, size(WT,2));

        % normalize so that we get a percentage difference like in the Figure 3
        % of Peng et al. (2012), http://dx.doi.org/10.1371/journal.pone.0034163
        % for example
        diff = WT - baselineVector;
        WTout = (diff ./ baselineVector);
        
        whos
        

    
    %% DEBUG
    if debugOn == 1
        
        scrsz = handles.style.scrsz;
        fig = figure('Color', 'w');
            set(fig, 'Position', [0.04*scrsz(3) 0.25*scrsz(4) 0.9*scrsz(3) 0.51*scrsz(4)])    
        
        % debug a point in 2D matrix
        fInd = 112;
        tInd = 188;
        
        whos
        
        input = WT(fInd, tInd)
        norm = baselineVector(fInd)
        difference = diff(fInd, tInd)
        output = WTout(fInd, tInd)
        outputPoint = 100 * (input - norm) / norm
        
        contourLevels = 256;
        % figure('color','w')
        ind = 1;
        sp(ind) = subplot(2,3,ind);
        contourf(WT, contourLevels, 'EdgeColor', 'none'); title('WT in')
        colorbar
        epochLimits(ind,:) = [min(min(WT)) max(max(WT))];
        
            hold on
            x = indicesTrim(1);
            y = 0;
            w = indicesTrim(2) - indicesTrim(1);
            h = length(f) - y;
            rr = rectangle('Position',[x y w h]);
            set(rr, 'LineWidth', 2, 'LineStyle', '--')
            tt = text(x+5, y+20, 'Baseline');
            
            hold off
        
        ind = ind + 1;
        sp(ind) = subplot(2,3,ind);
        contourf(WT_w_COI, contourLevels, 'EdgeColor', 'none'); title('WT in with COI')        
        colorbar        
        epochLimits(ind,:) = [min(min(WT_w_COI)) max(max(WT_w_COI))];
        
        ind = ind + 1;
        sp(ind) = subplot(2,3,ind);
        hold on
        plot(WT(fInd, :), 'Color', [0.2 0.6 1])            
        plot(round(mean([indicesTrim(1) indicesTrim(2)])), norm, 'o', 'MarkerFaceColor', 'k')            
        plot(diff(fInd, :), 'Color', [1 .2 .8])
        plot(WTout(fInd, :), 'Color', 'k')        
        
            leg = legend(['IN fInd=', num2str(fInd), ', f=', num2str(f(fInd)), ' Hz'], ['mean = ', num2str(norm)], ...
                         'Diff (WTin - baseline)', 'Out (diff ./ baseline');
            set(leg, 'FontSize', 7)
            legend('boxoff')
            title('Power of one frequency')
        
        ind = ind + 1;
        sp(ind) = subplot(2,3,ind);
        contourf(baselineVector, contourLevels, 'EdgeColor', 'none'); title('baselineVector')
        colorbar
        epochLimits(ind,:) = [min(min(baselineVector)) max(max(baselineVector))];
        
        ind = ind + 1;
        sp(ind) = subplot(2,3,ind);
        contourf(diff, contourLevels, 'EdgeColor', 'none'); title('Difference (WT-baseline)')
        colorbar
        epochLimits(ind,:) = [min(min(diff)) max(max(diff))];
        
        ind = ind + 1;
        sp(ind) = subplot(2,3,ind);
        contourf(WTout, contourLevels, 'EdgeColor', 'none'); title('Normalized Out (Diff ./ baseline)')
        epochLimits(ind,:) = [min(min(WTout)) max(max(WTout))];
        colorbar
        
        for i = 1 : ind-1
            caxis(sp(i), [min(min(epochLimits(1:ind-1,:))) max(max(epochLimits(1:ind-1,:)))]);
        end
        %set(sp(3), 'YLim', [min(min(epochLimits(1:ind-1,:))) max(max(epochLimits(1:ind-1,:)))])
        
      % AUTO-SAVE FIGURE
        try
            if handles.figureOut.debugON == 1 
                drawnow
                dateStr = plot_getDateString(); % get current date as string                
                fileNameOut = ['timeFreqNormalization_debug_', dateStr];
                disp(['         ... saving figure to disk (', fileNameOut, '.png]'])
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)                
            end
        catch err
            err
        end  
        
    end
    
    
    
    