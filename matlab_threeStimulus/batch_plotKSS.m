function batch_plotKSS(KSS, KSS_z, KSS_Norm, KSS_z_Norm, parameters, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempKSSPlot.mat';
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

    scrsz = get(0,'ScreenSize'); % get screen size for plotting
    fig = figure('Color','w','Name','KSS');
        set(fig, 'Position', [0.05*scrsz(3) 0.0625*scrsz(4) 0.92*scrsz(3) 0.9*scrsz(4)])  
    
        rows = 4;
        cols = 3;
        
        % sessions, after 2nd session, no KSS was asked
        t = [1 3 4];
        
        ind = 1;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'dark';
            [p, legStr] = batch_KSS_subplot(t, KSS.(condition), sp(ind), condition, ind, parameters, handles);
        
        ind = 2;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'dim';
            [p, legStr] = batch_KSS_subplot(t, KSS.(condition), sp(ind), condition, ind, parameters, handles);
        
        ind = 3;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'bright';
            [p, legStr] = batch_KSS_subplot(t, KSS.(condition), sp(ind), condition, ind, parameters, handles);
                    
            leg = legend(legStr);
                set(leg,'Position',[0.918996054333765 0.714130434782609 0.0443078913324709 0.199869565217391]);
                %legend('boxoff')
          
        ind = 4;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'dark';
            [p, legStr] = batch_KSS_subplot(t, KSS_z.(condition), sp(ind), condition, ind, parameters, handles);
        
        ind = 5;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'dim';
            [p, legStr] = batch_KSS_subplot(t, KSS_z.(condition), sp(ind), condition, ind, parameters, handles);
        
        ind = 6;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'bright';
            [p, legStr] = batch_KSS_subplot(t, KSS_z.(condition), sp(ind), condition, ind, parameters, handles);
            
        ind = 7;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'dark';
            [p, legStr] = batch_KSS_subplot(t, KSS_z_Norm.(condition), sp(ind), condition, ind, parameters, handles);
        
        
        ind = 8;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'dim';
            [p, legStr] = batch_KSS_subplot(t, KSS_z_Norm.(condition), sp(ind), condition, ind, parameters, handles);
        
        
        ind = 9;
        sp(ind) = subplot(rows,cols,ind);
        
            condition = 'bright';
            [p, legStr] = batch_KSS_subplot(t, KSS_z_Norm.(condition), sp(ind), condition, ind, parameters, handles);
        
        
                
        try
            if handles.figureOut.ON == 1    
                drawnow
                dateStr = plot_getDateString(); % get current date as string
                %cd(path.outputFigures)            
                fileNameOut = ['plot_KSS_', dateStr];
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)
                %cd(path.code)
            end
        catch err            
            if strcmp(err.identifier, '????')
                str = sprintf('%s\n%s', 'Crashing probably because you have not installed export_fig from Matlab File Exchange!', ...
                              'Download it from: http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig, and "File -> Set Path -> Add Folder"');
                error(str)
            else
                err
                err.identifier
            end
        end
        
    function [p, legStr] = batch_KSS_subplot(t, KSS, sp, condition, ind, parameters, handles)
        
        subjects = fieldnames(KSS);        
        ColorSet = distinguishable_colors(length(subjects));
            set(sp, 'ColorOrder', ColorSet);
                
        for sub = 1 : length(subjects)            
            y(sub,:) = KSS.(subjects{sub});
            legStr{sub} = subjects{sub};                    
        end
        
        p = plot(t, y, '-o');
        
        % style
        
            for sub = 1 : length(subjects)      
               set(p(sub), 'MarkerEdgeColor', 'none', 'MarkerFaceColor', ColorSet(sub,:), 'Color', ColorSet(sub,:))
            end
        
        
                
        lab(1) = xlabel('Session');
        lab(2) = ylabel('KSS');
        
        tit = title(condition);
        
        set(sp, 'XTick', t)
        
        set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1)
        set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold')
        set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')
        