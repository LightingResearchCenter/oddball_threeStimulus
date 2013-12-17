function batch_plotIntensityComparisonMAIN(statsOut, matricesSessionNorm, statsPer, erpComponent, erpFilterType, fieldValue, fileNameFields, stimulusType, chsToPlot, subjects, outlierOut, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempBatchPlot.mat';
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
    
    %{
    whos
    statsOut
    statsOut.darkCondition
    statsOut.darkCondition.target
    statsOut.darkCondition.target.dark
    statsOut.darkCondition.target.dark.Fz
    statsOut.darkCondition.target.dark.Fz.mean
    %}
    
    % Data "hidden" in structure        
    normalizationTypes = fieldnames(statsOut);
        erpTypes = fieldnames(statsOut.darkCondition);
            conditionTypes = fieldnames(statsOut.darkCondition.standard);
                chNames = fieldnames(statsOut.darkCondition.standard.dark);
                    statFields = fieldnames(statsOut.darkCondition.standard.dark.Cz);
                        noOfSessions = length(statsOut.darkCondition.standard.dark.Cz.mean);
    
    
    %% PLOT 1: Absolute / Normalization
        
        fig1 = figure('Color', 'w');
            set(fig1, 'Name', 'MEAN PLOT: Subjects averaged & Normalized')
            plot_componentAbsoluteNormalizedFigure(fig1, statsOut, normalizationTypes, erpComponent, erpFilterType, chsToPlot, erpTypes, fieldValue, subjects, handles)
            
    
    %% PLOT 2:        

        fig2 = figure('Color', 'w');
            chSelected = 1; % too clogged if you plot 2 channels to same figure
            plot_componentScatterFigure(fig2, statsOut, matricesSessionNorm, noOfSessions, normalizationTypes, erpComponent, erpFilterType, chsToPlot, chSelected,  erpTypes, fieldValue, subjects, handles)

        fig3 = figure('Color', 'w');
            chSelected = 2; % too clogged if you plot 2 channels to same figure
            plot_componentScatterFigure(fig3, statsOut, matricesSessionNorm, noOfSessions, normalizationTypes, erpComponent, erpFilterType, chsToPlot, chSelected, erpTypes, fieldValue, subjects, handles)

  