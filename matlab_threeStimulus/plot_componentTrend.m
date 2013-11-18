function plot_componentTrend(analyzed_target, analyzed_distracter, analyzed_standard, plotType, handles)

    %% DEBUG
    debugMatFileName = 'tempTrendPlot.mat';
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

    
    handles.style.fontGrey = [0.2 0.2 0.2];
    style = handles.style;
    
    
    %%     
    chNames = fieldnames(analyzed_target.ERPcomponents);
    %{
    analyzed_target
    analyzed_distracter
    analyzed_distracter.ERPcomponents
    analyzed_distracter.ERPcomponents
    analyzed_distracter.ERPcomponents.Fz
    analyzed_distracter.ERPcomponents.Fz{1}
    analyzed_distracter.ERPcomponents.Fz{1}.P3
    analyzed_distracter.ERPcomponents.Fz{1}.P3.meanAmplit
    %}
    
    %% CONDITION DATA
    
        if strcmp(plotType, 'allSamples')        
            consecutiveSamplesToAverage = 1;            
        elseif strcmp(plotType, 'averaged')
            consecutiveSamplesToAverage = handles.parameters.plotComponentAveragedSamples;            
        else
            error('Typo with the plotType variable?')
        end
        
        targetMAT = plot_convComponentCellsToMat(analyzed_target.ERPcomponents, consecutiveSamplesToAverage, handles);
        distracterMAT = plot_convComponentCellsToMat(analyzed_distracter.ERPcomponents, consecutiveSamplesToAverage, handles);
        standardMAT = plot_convComponentCellsToMat(analyzed_standard.ERPcomponents, consecutiveSamplesToAverage, handles);
      
    
    %%
    scrsz = handles.style.scrsz;
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.1*scrsz(3) 0.05*scrsz(4) 0.85*scrsz(3) 0.91*scrsz(4)])
    
        rows = 4;
        cols = 3;
    
        for i = 1 : handles.parameters.EEG.nrOfChannels
           
            % target
            ij = 1;
            ind = (i-1)*cols + ij;
            sp(i,ij) = subplot(rows,cols,ind);
                [p(i,ij,:), pNaN(i,ij), yLimits(i,ij,:)] = plot_trendComponents_perChannelSubplot(targetMAT, i, ij, 'target', chNames, rows, cols, consecutiveSamplesToAverage, handles);
            
            % distracter
            ij = 2;
            ind = (i-1)*cols + ij;
            sp(i,ij) = subplot(rows,cols,ind);
                [p(i,ij,:), pNaN(i,ij), yLimits(i,ij,:)] = plot_trendComponents_perChannelSubplot(distracterMAT, i, ij, 'distracter', chNames, rows, cols, consecutiveSamplesToAverage, handles);
            
            % standard
            ij = 3;
            ind = (i-1)*cols + ij;
            sp(i,ij) = subplot(rows,cols,ind);
                [p(i,ij,:), pNaN(i,ij), yLimits(i,ij,:)] = plot_trendComponents_perChannelSubplot(standardMAT, i, ij, 'standard', chNames, rows, cols, consecutiveSamplesToAverage, handles);
            
        end
        
        yLimitMin = min(min(min(yLimits)));
        yLimitMax = max(max(max(yLimits)));
        set(sp, 'YLim', [yLimitMin yLimitMax])
        set(sp, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-2) % , 'Color', handles.style.fontGrey)
    
        
    
    
    %% SUBFUNCTIONS
    
    function [p, pNaN, yLimits] = plot_trendComponents_perChannelSubplot(matrixIn, i, ij, stimulusType, chNames, rows, cols, consecutiveSamplesToAverage, handles)

        style = handles.style;
        parameters = handles.parameters;
        
        %% GET INPUT DATA
        
            componentNames = fieldnames(matrixIn.(chNames{i}));        
            if ~strcmp(componentNames, 'RT')
                statFields = fieldnames(cellIn.(chNames{ch}){epochs}.(componentNames{component}));
            else
                statFields = {'only RT here, thus length is 1'};
            end

            % easier variable names            
            N2 = matrixIn.(chNames{i}).N2.peakMeanAmplit.vector;
            P3 = matrixIn.(chNames{i}).P3.peakMeanAmplit.vector;
            P3_N2 = matrixIn.(chNames{i}).P3_N2.peakMeanAmplit.vector;
            
            if consecutiveSamplesToAverage ~= 1
                N2_SD = matrixIn.(chNames{i}).N2.peakMeanAmplit.sdOut;
                P3_SD = matrixIn.(chNames{i}).P3.peakMeanAmplit.sdOut;
                P3_N2_SD = matrixIn.(chNames{i}).P3_N2.peakMeanAmplit.sdOut;
            end
            
            n = sum(~isnan(P3_N2));
                % check later if this actually is correct, and when do the
                % outlier/artifact rejection occur, correct to
                % analyze_getERPcomponents() if needed
            
            nanIndices = isnan(P3_N2);
                
            x = linspace(1,length(N2),length(N2))';

        %% PLOT
        
            hold on
            if consecutiveSamplesToAverage == 1
                p(1:3) = plot(x, N2, 'o', x, P3, 'o', x, P3_N2, 'o');
            else
                p(1) = errorbar(x,N2,N2_SD, 'ro');
                p(2) = errorbar(x,P3,P3_SD, 'go');
                p(3) = errorbar(x,P3_N2,P3_N2_SD, 'bo');
            end
            
            
            % Style
            
                if i == rows
                    xLab = xlabel('Trial #');
                else
                    xLab = xlabel(' ');
                end
                
                if ij == 1
                    yLab = ylabel('\muV');
                else
                    yLab = ylabel(' ');
                end
        
                yLimits = [min(min([N2 P3 P3_N2])) max(max([N2 P3 P3_N2]))];
                
                    % PLOT the NaNs
                    yNaN = repmat(yLimits(1), sum(nanIndices == 1), 1);
                    xNaN = x(nanIndices);
                    if isempty(xNaN) || isempty(yNaN)
                        pNaN = plot(NaN, NaN); % no values when averaging trials typically
                    else
                        pNaN = plot(xNaN, yNaN, 'ko');
                    end
                    hold off
                    
                                        
                if ij == 1 && i == 1
                    leg = legend('N2', 'P3', 'P3-N2', 'NaN');
                        set(leg, 'Position',[0.323652694610778 0.869775541795667 0.104790419161677 0.109713622291022])
                        legend('boxoff')    
                else
                    leg = legend('');
                        legend('hide')
                
                end
                
                if consecutiveSamplesToAverage == 1
                    style.markerSize = style.markerSize - 2;
                end
                
                set(p, 'MarkerSize', style.markerSize, 'MarkerEdgeColor', [.4 .4 .4], 'Color', style.fontGrey)
                    set(p(1), 'MarkerFaceColor', 'r')
                    set(p(2), 'MarkerFaceColor', 'g')
                    set(p(3), 'MarkerFaceColor', 'b')

                set(pNaN, 'MarkerSize', style.markerSize, 'MarkerEdgeColor', [.4 .4 .4])
                    set(pNaN, 'MarkerFaceColor', 'k')
                
                if i == 1
                    
                    if strcmp(stimulusType, 'standard')
                        multiplier = parameters.oddballTask.nrOfStandardsPerCycle;
                        freq = parameters.oddballTask.standardFreq;
                    elseif strcmp(stimulusType, 'target')
                        multiplier = 1;
                        freq = parameters.oddballTask.targetFreq;
                    elseif strcmp(stimulusType, 'distracter')
                        multiplier = 1;
                        freq = parameters.oddballTask.distractFreq;
                    else
                        multiplier = 1;
                    end
                    
                    titStr = sprintf('%s\n%s\n%s', [stimulusType, ' (', parameters.BioSemi.chName{i + parameters.BioSemi.chOffset}, ')'], ...
                                    ['freq: ', num2str(freq), ' Hz'], ...
                                    ['n = ', num2str(n), ' / ', num2str(multiplier*parameters.oddballTask.repeatsOfCycle / consecutiveSamplesToAverage)]);
                    tit = title(titStr);
                else
                    tit = title(['(', parameters.BioSemi.chName{i + parameters.BioSemi.chOffset}, ')']);
                end             
                
                set([xLab yLab], 'FontName', style.fontName, 'FontSize', style.fontSizeBase-2, 'FontWeight', 'bold', 'Color', handles.style.fontGrey)
                set(tit, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-1, 'FontWeight', 'bold', 'Color', handles.style.fontGrey)
                set(leg, 'FontName', style.fontName, 'FontSize', style.fontSizeBase-1, 'FontWeight', 'bold', 'Color', handles.style.fontGrey)
                
    
    function matrixOut = plot_convComponentCellsToMat(cellIn, consecutiveSamplesToAverage, handles)
        
        chNames = fieldnames(cellIn);
        
        for ch = 1 : length(chNames)
           
            for epochs = 1 : length(cellIn.(chNames{ch}))
                
                % cellIn.(chNames{ch}){epochs}
                componentNames = fieldnames(cellIn.(chNames{ch}){epochs});
                
                for component = 1 : length(componentNames)
                   
                    % cellIn.(chNames{ch}){epochs}.(componentNames{component})
                    % componentNames
                    % [ch epochs component]
                    
                    if ~strcmp(componentNames{component}, 'RT')
                        statFields = fieldnames(cellIn.(chNames{ch}){epochs}.(componentNames{component}));
                    else
                        statFields = {'only RT here, thus length is 1'};
                    end
                    
                        for statParam = 1 : length(statFields)                            
                            if ~strcmp(componentNames{component}, 'RT')
                                matrixOut.(chNames{ch}).(componentNames{component}).(statFields{statParam}).vector(epochs) = cellIn.(chNames{ch}){epochs}.(componentNames{component}).(statFields{statParam});
                            else
                                matrixOut.(chNames{ch}).(componentNames{component}).vector(epochs) = cellIn.(chNames{ch}){epochs}.(componentNames{component});
                            end
                                
                        end
                    
                end
                
            end
            
        end
        
        if consecutiveSamplesToAverage == 1           
           % nothing
            
        else
            for ch = 1 : length(chNames)
                for component = 1 : length(componentNames)
                    if ~strcmp(componentNames{component}, 'RT')
                        statFields = fieldnames(cellIn.(chNames{ch}){epochs}.(componentNames{component}));
                    else
                        statFields = {'only RT here, thus length is 1'};
                    end
                    
                    for statParam = 1 : length(statFields)
                        if ~strcmp(componentNames{component}, 'RT')
                            
                            ik = 1;
                            averageCounter = 1;
                            
                            for ij = 1 : length(matrixOut.(chNames{ch}).(componentNames{component}).(statFields{statParam}).vector)                                                             
                                tempVector(ik) = matrixOut.(chNames{ch}).(componentNames{component}).(statFields{statParam}).vector(ij);
                                ik = ik + 1;
                                if rem(ij,consecutiveSamplesToAverage) == 0
                                    % ik
                                    % tempVector
                                    meanValue(averageCounter) = nanmean(tempVector);
                                    sdValue(averageCounter) = nanstd(tempVector);
                                    nValue(averageCounter) = sum(~isnan(tempVector));
                                    ik = 1; % reset
                                    averageCounter = averageCounter + 1;
                                end
                            end                       
                            
                            % a = 'components'
                            % meanValue
                            % sdValue
                            matrixOut.(chNames{ch}).(componentNames{component}).(statFields{statParam}).vector = meanValue;
                            matrixOut.(chNames{ch}).(componentNames{component}).(statFields{statParam}).sdOut = sdValue;
                            matrixOut.(chNames{ch}).(componentNames{component}).(statFields{statParam}).n = nValue;
                            
                        else
                            
                            ik = 1;
                            averageCounter = 1;
                            
                            for ij = 1 : length(matrixOut.(chNames{ch}).(componentNames{component}).vector)                                                             
                                tempVector(ik) = matrixOut.(chNames{ch}).(componentNames{component}).vector(ij);
                                ik = ik + 1;
                                if rem(ij,consecutiveSamplesToAverage) == 0
                                    meanValue(averageCounter) = nanmean(tempVector);
                                    sdValue(averageCounter) = nanstd(tempVector);
                                    nValue(averageCounter) = sum(~isnan(tempVector));
                                    ik = 1; % reset
                                    averageCounter = averageCounter + 1;
                                end
                            end
                            
                            a = 'RT'
                            meanValue
                            sdValue
                            matrixOut.(chNames{ch}).(componentNames{component}).vector = meanValue;
                            matrixOut.(chNames{ch}).(componentNames{component}).sdOut = sdValue;
                            matrixOut.(chNames{ch}).(componentNames{component}).n = nValue;
                            
                        end

                    end
                end
            end
            
        end
                    
        
        