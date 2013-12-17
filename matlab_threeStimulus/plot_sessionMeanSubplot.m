function [p, styleHandles, yLims] = ...
    plot_sessionMeanSubplot(statData, chToPlot, normalizationTypes, erpTypes, rows, cols, normType, stim, ch, index, fieldValue, erpComponent, erpDataType, subjects, handles)

    i = 1;

    x = 1:1:4;

    % displace the x-axis a bit for better visualization
    x_displ = 0.175;
    x1 = x;
    x2 = x - x_displ;
    x3 = x + x_displ;

    % statData
    % statData.target
    % erpTypes{stim}
    % statData.(erpTypes{stim})
    if strcmp(erpComponent, 'RT')       
        if strcmp(normalizationTypes{normType}, 'absolute')
            secondsToMilliseconds = 1000;
        else
            secondsToMilliseconds = 1;
        end            
        try
            y1 = secondsToMilliseconds * statData.(erpTypes{stim}).dark.(chToPlot).mean;
        catch err
            err
            statData
            erpTypes{stim}
            statData.(erpTypes{stim})
            statData.(erpTypes{stim}).dark
            statData.(erpTypes{stim}).dark.(chToPlot)
            error('indexing incorrect')
        end
        y2 = secondsToMilliseconds * statData.(erpTypes{stim}).dim.(chToPlot).mean;
        y3 = secondsToMilliseconds * statData.(erpTypes{stim}).bright.(chToPlot).mean;

        % these errors are SDs of the means, thus for normalized values
        % you won't have by definition any error, add later the LE and
        % UE obtained in batch_normalizeComponents() using the
        % low-level subfunction batch_normalizeLowLevel()
        err1 = secondsToMilliseconds * statData.(erpTypes{stim}).dark.(chToPlot).SD;
        err2 = secondsToMilliseconds * statData.(erpTypes{stim}).dim.(chToPlot).SD;
        err3 = secondsToMilliseconds * statData.(erpTypes{stim}).bright.(chToPlot).SD;
    else
        % plotStatsDebug = statData.(erpTypes{stim}).dark.(chToPlot)
        y1 = statData.(erpTypes{stim}).dark.(chToPlot).mean;
        y2 = statData.(erpTypes{stim}).dim.(chToPlot).mean;
        y3 = statData.(erpTypes{stim}).bright.(chToPlot).mean;

        err1 = statData.(erpTypes{stim}).dark.(chToPlot).SD;
        err2 = statData.(erpTypes{stim}).dim.(chToPlot).SD;
        err3 = statData.(erpTypes{stim}).bright.(chToPlot).SD;
    end

    % PLOT
    hold on
    p(1) = errorbar(x1, y1, err1, 'ko');
    p(2) = errorbar(x2, y2, err2, 'r*');
    p(3) = errorbar(x3, y3, err3, 'ro');
    hold off

    % Annotations

        % number of non-NaN subjects per plot
        %{
        yOffset = 1.1;
        yOffset2 = 0.01;
        for ll = 1 : length(x1)
            noOfNotNans = sum(~isnan(y1(ll,:)));
            tx(1,ll) = text(x1(ll), max(y1(ll,:))*yOffset + yOffset2, num2str(noOfNotNans));
        end
           
        for ll = 1 : length(x2)
            noOfNotNans = sum(~isnan(y3(ll,:)));
            tx(2,ll) = text(x2(ll), max(y2(ll,:))*yOffset + yOffset2, num2str(noOfNotNans));
        end
        
        for ll = 1 : length(x3)
            noOfNotNans = sum(~isnan(y2(ll,:)));
            tx(3,ll) = text(x3(ll), max(y3(ll,:))*yOffset + yOffset2, num2str(noOfNotNans));
            
        end
        set(tx, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontAngle', 'Italic')
        %}
    
        % title & ylabel
        if rem(index-1,cols) == 0 || index == 1 % for each row
            normString = normalizationTypes{normType};

            if strcmp('absolute', normString)
                if strcmp(erpComponent, 'RT')
                    yLabelString = 'Latency [ms]';
                else
                    yLabelString = '\muV';
                end
            else
                yLabelString = '\Delta';
            end                
        elseif (index - (floor(index/cols)*cols)) == 2
            normString = erpComponent;
            yLabelString = ' ';
        elseif (index - (floor(index/cols)*cols)) == 4
            normString = fieldValue;                
            yLabelString = ' ';
        elseif (index - (floor(index/cols)*cols)) == 0
            normString = erpDataType;                
            yLabelString = ' ';                
        else                
            normString = ' ';
            yLabelString = ' ';
        end

        if rem(index, length(chToPlot)) ~= 0
            stimString = erpTypes{stim};
        else
            stimString = ' ';
        end
        if strcmp(erpComponent, 'RT')
            titString = sprintf('%s', normString);    
        else
            titString = sprintf('%s\n%s\n%s', normString, stimString, chToPlot);    
        end
        styleHandles.tit = title(titString);

        % xlabel
        if (rows - 1)*cols < index
            xLabString = sprintf('%s', 'Session');
        else
            xLabString = ' ';
        end
        styleHandles.xLab = xlabel(xLabString);

        % ylabel
        styleHandles.yLab = ylabel(yLabelString);

        % LEGEND
        if strcmp(erpComponent, 'RT')
            if index == 1
                leg = legend('dark', '10lux', '40lux', 3);
                    %legend('boxoff')
                    set(leg, 'Location', 'NorthWest')
                    set(leg, 'Position', [0.807142193308545 0.911998913043478 0.160780669144981 0.0654891304347826])
            else
                leg = legend(' ', ' ', ' ', 3);
                    legend('hide')
            end

        else
            if (index - (floor(index/cols)*cols)) == 0
                leg = legend('dark', '10lux', '40lux', 3);
                    %legend('boxoff')
                    set(leg, 'Location', 'NorthWest')
                    set(leg, 'Position', [0.9034 0.2538 0.050776 0.04918])
            else
                leg = legend(' ', ' ', ' ', 3);
                    legend('hide')
            end
        end


    % Style
    set(gca, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)
    set(styleHandles.tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
        set(styleHandles.tit, 'HorizontalAlignment', 'center')
    set(styleHandles.xLab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
    set(styleHandles.yLab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2)
        set(leg, 'EdgeColor', [.4 .4 .4])

    % plots
    set(p(1),'MarkerFaceColor',[0 0 0],'Color',[0 0 0], 'Marker', 's');
    set(p(2),'MarkerFaceColor',[0 0.20 0.80],...
        'Marker','pentagram',...
        'Color',[0 0.20 0.80]);
    set(p(3),'MarkerFaceColor',[1 0 0],'Color',[1 0 0]);
    set(p, 'MarkerSize', handles.style.markerSize+1)

    % get y limits
    yLims = [min(min(([y1 y2 y3]))) max(max(([y1 y2 y3])))];

