function [dataMatrixOut, outlierIndices] = batch_excludeOutliersDuringBatch(dataMatrix, subjects, callFromWhere, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempBatchExclusion.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
            debugPlot = 1;
        else
            if handles.flags.saveDebugMATs == 1
                path = handles.path;
                save('debugPath.mat', 'path')
                save(fullfile(path.debugMATs, debugMatFileName))            
            end
            debugPlot = 0;
        end
    end
            
    
    if strcmp(callFromWhere, 'normalize')
        
    elseif strcmp(callFromWhere, 'average')
       
    else
        callFromWhere
    end
    
    subjects = unique(subjects);
    [noOfEpochsOrDataPoints, noOfSubjects] = size(dataMatrix);
        
    % calculate parameters of the input
    x_mean = nanmean(dataMatrix);
    x_median = nanmedian(dataMatrix);
    x_sd = nanstd(dataMatrix);
    UE = x_mean + x_sd;
    LE = x_mean - x_sd;
    
    % outlier rejection
    sigmaMultiplier = 2;
    outlierExcl_UE = x_mean + (sigmaMultiplier * x_sd);
    outlierExcl_LE = x_mean - (sigmaMultiplier * x_sd);
    
    % reject
    dataMatrixOut = dataMatrix;
    dataMatrixOut(dataMatrix > outlierExcl_UE) = NaN;
    dataMatrixOut(dataMatrix < outlierExcl_LE) = NaN;
    outlierIndices = isnan(dataMatrixOut);
    
    if debugPlot == 1
        scrsz = get(0,'ScreenSize'); % get screen size for plotting
        fig = figure('Color', 'w');                

            if strcmp(callFromWhere, 'average')
                
                set(fig, 'Position', [0.06*scrsz(3) 0.34*scrsz(4) 0.9*scrsz(3) 0.6*scrsz(4)])    
                rows = 2;
                cols = 3;
                
                    ind = 1;
                    sp(ind) = subplot(rows,cols,ind);
                        p1 = plot(dataMatrix);
                        lab(ind,1) = xlabel('Epoch no'); lab(ind,2) = ylabel('Abs. Value'); tit(ind) = title('All Subjects / All Epochs');

                    ind = 2;
                    sp(ind) = subplot(rows,cols,ind);
                        p2 = plot(nanmean(dataMatrix), 'ko');
                        lab(ind,1) = xlabel('Subject'); lab(ind,2) = ylabel('Mean Value'); tit(ind) = title('All Subjects (epochs averaged)');
                        set(p2, 'MarkerFaceColor', [0 .6 1])

                    ind = 3;
                    sp(ind) = subplot(rows,cols,ind);
                        p3 = plot(nanmedian(dataMatrix), 'ko');
                        lab(ind,1) = xlabel('Subject'); lab(ind,2) = ylabel('Mean Value'); tit(ind) = title('All Subjects (epoch median)');
                        set(p3, 'MarkerFaceColor', [1 .2 .6])
                    
            elseif strcmp(callFromWhere, 'normalize')
                
                set(fig, 'Position', [0.06*scrsz(3) 0.64*scrsz(4) 0.9*scrsz(3) 0.3*scrsz(4)])  
                rows = 1;
                cols = 3;
                
                    ind = 1;
                    sp(ind) = subplot(rows,cols,ind);
                        hold on
                        offs = 0.5;
                        p1 = plot(dataMatrix, 'ro');
                        xLims = get(gca, 'XLim');
                        
                        p2 = plot(mean(xLims)+offs, x_mean, 'ko');
                        p3 = plot(mean(xLims)-offs, x_median, 'ko');
                        p4(1) = line(xLims, [UE UE]);
                        p4(2) = line(xLims, [LE LE]);                        
                        
                        p5(1) = line(xLims, [outlierExcl_UE outlierExcl_UE]);
                        p5(2) = line(xLims, [outlierExcl_LE outlierExcl_LE]);
                        
                        p6 = plot(find(outlierIndices == 1), dataMatrix(outlierIndices), 'kx');
                        hold off
                        
                        lab(ind,1) = xlabel('Subject'); lab(ind,2) = ylabel('Norm. Value'); tit(ind) = title('All Subjects / All Epochs');
                        set(p1, 'MarkerFaceColor', [1 .6 1])
                        set(p2, 'MarkerFaceColor', [0 .6 1])
                        set(p3, 'MarkerFaceColor', [1 .2 .6])
                        set([p2 p3], 'MarkerSize', 7)
                        set(p6, 'MarkerSize', 13)
                        set(p4, 'Color', [.4 .4 .42], 'LineStyle', '--')
                        set(p5, 'Color', 'k', 'LineStyle', ':')
                        
                        leg = legend([p1 p2 p3 p4(1) p5(1)], 'Subjects', 'Mean', 'Median', 'SD', ['Outlier Thr (', num2str(sigmaMultiplier), '\sigma)']);
                            set(leg, 'Position',[0.0411706349206348 0.669444444444445 0.0677910052910053 0.275793650793651])
                                legend('boxoff')
                            
                        % annotate with subjects
                        yOffs = 0.1;
                        for i = 1 : noOfSubjects
                            tx(i) = text(i, dataMatrix(i), subjects{i}, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');
                        end
                        


            end
            
            set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
            set(sp, 'XLim', [1 max(xLims)]) 
            set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
            set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')  

        % sdd

    end