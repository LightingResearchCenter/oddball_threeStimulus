function plot_allTheEpochsToSingleSubplot(spHandle, t, EEGepochs, parameters, yOffset)

    numberOfEpochs = size(EEGepochs,3);
    xOff = 20; % [ms]

    hold on
    for ij = 1 : numberOfEpochs

       % axes(spHandle) 

       yOff = yOffset*ij; % update the horizontal line (baseline)

       y(ij,:,:) = EEGepochs(:,:,ij);

       % Horizontal line (baseline)
       l(ij) = line([min(t) max(t)], [yOff yOff], 'Color', 'k');                                     

       % Filtered
       p_ERP(ij, :) = plot(t, yOff + squeeze(y(ij,:,:)));

       yTickPos(ij) = yOff; % save for yTick positions (Trial)

       yTickLabel{ij} = num2str(ij);
       drawnow

       maxY(ij) = max(max(yOff + abs(squeeze(y(ij,:,:)))));
       maxAmplitude(ij) = maxY(ij) - yOff;
       minY(ij) = min(min(yOff + squeeze(y(ij,:,:))));

       maxHandle(ij) = text(max(t)+xOff, yOff, num2str(maxAmplitude(ij), '%3.0f'));
       if ij == numberOfEpochs
           yOff = yOffset * (ij + 1);
           maxHandle(ij+1) = text(max(t)+xOff, yOff, 'max');
       end

    end

    set(maxHandle, 'FontSize', 6')
    set(gca, 'YTick', yTickPos, 'YTickLabel', yTickLabel)
    yLims = get(gca, 'YLim');        
    yMin = nanmin(minY);
    yMax = nanmax(maxY);

    % not optimal for incoming data as gets the min and max from whole
    % matrix, fix later
    try
        set(gca, 'YLim', [yMin yMax])
    catch err            
        if strcmp(err.identifier, 'MATLAB:hg:propswch:PropertyError')
            warning('yLimits are the same, most likely all the epochs are considered as artifacts?')
        end
    end