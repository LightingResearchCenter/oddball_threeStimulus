function [sp_i, sp] = plot_CRAPandFIXED_steps(fig, sp, sp_i, leg, rows, cols, ...
                        NaN_indices_moving, NaN_indices_step, EEGfixedIndices, EOGfixedIndices, vDiffOutMovWindow, vDiffOutStep, ...
                        subplotIndices, parameters, handles)

    % whos
                    
    sp_i = sp_i + 1;
    i = 1;
    sp(sp_i) = subplot(rows,cols, subplotIndices(1));   
       
        % weigh differently for visualization
        offs = 0.2;
        for ch = 1 : size(EEGfixedIndices,2)
            EEGfixedIndices(:,ch) = (1 - ((ch-1)*offs)) * EEGfixedIndices(:,ch);
        end
    
        hold on
        p1a = plot(EEGfixedIndices, 'o','MarkerFaceColor',[.49 .49 .49]);
        p1b = plot(0.2*EOGfixedIndices, 'mo','MarkerFaceColor',[.49 .49 .49]);
        hold off        
    
        lab(i,1) = xlabel('Epochs');
        lab(i,2) = ylabel(['Artifact [~= 0]']);
        numberOfFixedArtifacts = sum(logical(EEGfixedIndices + EOGfixedIndices), 1);
        tit(i) = title(['Fixed Indices: [', num2str(numberOfFixedArtifacts), ']']);
        leg(1) = legend([p1a(1) p1b(1)], ['EEG chs, ', num2str(parameters.artifacts.fixedThr), ' \muV'], ['EOG ', num2str(parameters.artifacts.fixedThrEOG), ' \muV'], 'Location', 'Best');
            legend('boxoff')

    sp_i = sp_i + 1;
    i = 2;
    sp(sp_i) = subplot(rows,cols, subplotIndices(2));
    
        hold on        
        yThr = handles.parameters.artifacts.CRAP.movWind_ampTh(2); % Threshold
        p2 = plot(1:length(vDiffOutMovWindow), vDiffOutMovWindow, 1:length(vDiffOutMovWindow), yThr:yThr);   
        hold off
    
        lab(i,1) = xlabel('Epochs');
        lab(i,2) = ylabel(['max(diff) [\muV]']);
        tit(i) = title(['CRAP - MovWindow: [', num2str(sum(NaN_indices_moving,1)), ']']);
        
            leg(2) = legend(p2, 'Cz', 'Fz', 'Pz', 'Oz', 'Location', 'Best');
            legend('boxoff')

    sp_i = sp_i + 1;
    i = 3;
    sp(sp_i) = subplot(rows,cols, subplotIndices(3));
    
        yThr = handles.parameters.artifacts.CRAP.step_ampTh; % Threshold
        p3 = plot(1:length(vDiffOutStep), vDiffOutStep, 1:length(vDiffOutStep), yThr:yThr);        
        
    
        lab(i,1) = xlabel('Epochs');
        lab(i,2) = ylabel(['max(diff) [\muV]']);
        tit(i) = title(['CRAP - Step: [', num2str(sum(sum(NaN_indices_step,1))), ']']);
        
            leg(3) = legend('EOG', 'Location', 'Best');
            legend('boxoff')

    sp_i = sp_i + 1;
    i = 4;
    sp(sp_i) = subplot(rows,cols, subplotIndices(4));       
    
        yMov = sum(NaN_indices_moving,2)/parameters.EEG.nrOfChannels;
        yStep = logical(sum(NaN_indices_step,2));
        hold on
        b(1) = bar(yMov, 0.5, 'g', 'EdgeColor', 'none');
        b(2) = bar(yStep, 0.5, 'b', 'EdgeColor', 'none');
        hold off
        
            alpha = .55;
            ch = get(b(1),'child');
                set(ch,'facea',alpha)                
            ch = get(b(2),'child');
                set(ch,'facea',alpha)
        
        lab(i,1) = xlabel('Epochs');
        lab(i,2) = ylabel(['Artifact']);
        tit(i) = title('CRAP Artifacts');
        xlim([1 length(yMov)])
        ylim([0 1.2])
                
        leg(4) = legend(['MovWin, n=', num2str(sum(yMov ), '%3.2f')],...
                        ['Step, n=', num2str(sum(yStep), '%3.2f')], 'Location', 'Best');
        set(leg(4), 'Position', [0.443420186103556 0.274500986618277 0.104949874686717 0.0354143337066069])
        legend('boxoff')

    set(leg, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1) 
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold') 
    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2, 'FontWeight', 'bold') 
