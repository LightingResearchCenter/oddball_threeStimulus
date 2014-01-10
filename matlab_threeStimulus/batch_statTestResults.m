function stat_results = batch_statTestResults(inputMatrix, stat_assumptions, normFieldName, stimType, chName, conditions, sessions, subjects, parameters, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatTesting.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
            localCall = 1;
        else
            if handles.flags.saveDebugMATs == 1
                path = handles.path;
                save('debugPath.mat', 'path')
                save(fullfile(path.debugMATs, debugMatFileName))            
            end
            localCall = 0;
        end
    end

    parameters.stats.anova_pThreshold = 0.05;
    
    [noOfConditions, noOfSessions, noOfSubjects] = size(inputMatrix);
    
    condVect = conditions
    sesVect = sessions;
    subVect = unique(subjects);            
    
    %{
    normFieldName
    if strcmp(normFieldName, 'firstSession')
        goesHere = 1
        chName
        if strcmp(chName, 'Pz')
            goesHere = 2
            stimType
            if strcmp(stimType, 'target')
                fasa
            end
        end
    end
    %}
    
    
    %% Debug plot the input
    
        if localCall == 1
            scrsz = get(0,'ScreenSize'); % get screen size for plotting
            fig = figure('Color','w','Name', 'inputMatrix Scatter');
                set(fig, 'Position', [0.01*scrsz(3) 0.64*scrsz(4) 0.3*scrsz(3) 0.3*scrsz(4)])  
                rows = 3; cols = 1;
                xOffs = 0.15;
                x1 = (1:4) - xOffs;

            hold on
            for cond = 1 : noOfConditions           

                x = x1 + ((cond - 1) * xOffs);
                p(cond,:) = plot(x, squeeze(inputMatrix(cond,:,:)), 'o');
                titStr = sprintf('%s\n%s\n%s', normFieldName, [stimType, ', ', chName], '1st column of session = dark, then dim, and then bright');
                title(titStr)
                xlabel('Colors correspond to different subjects')

            end
            hold off
            figure
        end
    
    %% Vectorize the input 3D matrix
        
        % concenate the sessions and subjects together (3D -> 2D)
        conditionsLeft = reshape(inputMatrix, noOfConditions, noOfSessions*noOfSubjects);        

            % replicate the subjects group
            nk = size(subVect,2);
            for i = 1 : noOfSessions
                idx(i, :) = [1:nk];
            end
            subVect = subVect(:,idx(:));

            % replicate the session group
            nk = size(sesVect,1);
            idx2 = [1:nk]';        
            idx2 = repmat(idx2, noOfSubjects, 1);   
            sesVect = sesVect(idx2(:),:);        

        % Vectorize INPUT DATA (2D -> 1D)

            inputVector = conditionsLeft(:);

        % Replicate again the groups

            % replicate the subjects group        
            nk1 = size(subVect,2);
            idx3 = [1:nk1]';
            idx3 = repmat(idx3, noOfConditions, 1);
            subVect = (subVect(:, idx3(:)))';

            % replicate the session group        
            nk2 = size(sesVect,1);
            idx4 = [1:nk2]';        
            idx4 = repmat(idx4, noOfConditions, 1);  
            sesVect = sesVect(idx4(:),:);

            % replicate the condition group        
            nk3 = size(condVect,1);
            for i = 1 : noOfSubjects*noOfSessions
                idx5(i,:) = [1:nk3];
            end
            condVect = condVect(idx5(:),:);

        % Combine the groups

            group = {condVect sesVect};
            %group = {condVect sesVect subVect};
    

    %% If Data are normal (Gaussian) and deviations are homogeneous then you can do an ANOVA 
    
        % --> significant if p < 0.05

        % Fix this if to actually take into an account what the assumption
        % test gave you about whether you can do ANOVA or not
        
        if 1 == 1 % sw.W(1) < 1 && sw.W(2) < 1 && sw.W(2) < 3 && btestOutGrouped{i}.P > handles.bartlett_pThreshold                

            % The ANOVA test makes the following assumptions about the data in X:
                % All sample populations are normally distributed.
                % All sample populations have equal variance.
                % All observations are mutually independent.                

            % the Matlab Statistics Toolbox                
            alpha  = parameters.stats.anova_pThreshold;
            
            if localCall == 1
                displayMode = 'on';
            else
                displayMode = 'off';
            end
            
            % For an introduction to n-Anova see e.g.
            % http://www.mathworks.com/help/stats/anova.html#bqttd20-1
            
            [p,table,stats,terms] = anovan(inputVector, group, 'model','interaction', ...
                                    'varnames',{'Condition';'Session'}, 'display', displayMode);

                [c,m,h,nms] = multcompare(stats,'display', displayMode);
                [nms num2cell(m)]
                % Notches in the boxplot provide a test of group medians 
                % (see boxplot) different from the F test for means in
                % the ANOVA table             
                
                % package to output
                stat_results.anovan.p = p;
                stat_results.anovan.table = table;
                stat_results.anovan.stats = stats;
                stat_results.anovan.terms = terms;
           
        else
            
            disp('   ASSUMPTIONS ARE NOT TRUE so that you could do ANOVA')
                
        end