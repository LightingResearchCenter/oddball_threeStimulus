function statsAdv = stat_signifTestWrapper(inputMatrix, parameters, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatSignifWrapper.mat';
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
       
   n = length(inputMatrix);

   if handles.displayMidResults == 1
       disp(' ');        
       disp('    ======= '); 
       disp(' ')
   end

    %% 1) You should test for the normality of distribution of the data. (Shapiro-Wilk normality test) 
    % ---> ok if p > 0.05

        % IMPLEMENTATION by: Ahmed Ben Saïda
        % http://www.mathworks.com/matlabcentral/fileexchange/13964
        % Shapiro-Wilk parametric hypothesis test of composite normality, for sample size 3<= n <= 5000. 
        % Based on Royston R94 algorithm. 
        % This test also performs the Shapiro-Francia normality test for platykurtic samples.
        cd(handles.path.statFunc)

        % inputs
        X = X_grouped{4}.paired(:,i)
        alpha = handles.shapWilk_pThreshold;
        tail = 1; % default value, check swtest.m for explations

        % call the subfunction
        if n >= 3
            [sw.H(i), sw.pValue(i), sw.W(i)] = swtest(X, alpha, tail);
        else
            sw.H(i) = NaN;
            sw.pValue(i) = NaN;
            sw.W(i) = NaN;
        end

        if handles.displayMidResults == 1
            disp(['     ', 'Shapiro-Wilk: [swtest.m]'])
            disp(['         H :', num2str(sw.H(i))])
            disp(['         p :', num2str(sw.pValue(i))])
            disp(['         W :', num2str(sw.W(i))])
        end

            % NULL hypothesis is that the distribution is NOT Gaussian so
            % if the p-value is higher than your alpha, then the
            % distribution is Gaussian, i.e. if W < 1

        cd(currDir)

    %% 2) If valid you test homogeneity of distribution of the deviation. (Bartlett's K-squared test) 
    % --> ok if p > 0.05

        if sw.W(i) < 1 % if the distribution is indeed GAUSSIAN

            cd(handles.path.statFunc)

            % Use the subfunction provided by Antonio
            % Trujillo-Ortiz to reduce dependency on additional
            % toolboxes such as Statistics Toolbox that has the
            % function "barttest" (?)
            % http://www.mathworks.com/matlabcentral/fileexchange/3314-btest                                       
            data   = X_grouped{4}.paired(:,i);
            sample = ones(length(data),1);
            X      = [data sample];
            alpha  = handles.bartlett_pThreshold;

            % first with Trujillo-Ortiz's subfunction
            btestOut(i) = Btest(X,alpha);

        if handles.displayMidResults == 1
            disp(['     ', 'Bartlett K-squared test: [btest.m]'])                    
            disp(['         variance      :', num2str(btestOut(i).var)])
            disp(['         deg. freedom  :', num2str(btestOut(i).v)])
            disp(['         p-probability :', num2str(btestOut(i).P)])
            disp(['         chi^2         :', num2str(btestOut(i).X2)])
            disp(['         F-statistic   :', num2str(btestOut(i).F)]); 
            disp(' ');
        end
            cd(currDir)

        else
            if handles.displayMidResults == 1
                disp(['     WARNING: The distribution WAS NOT Gaussian for the "', origGroup{i}.header, '!"']); disp(' ');
            end
            btestOut(i).var = NaN;
            btestOut(i).v   = NaN;
            btestOut(i).P   = NaN;
            btestOut(i).X2  = NaN;
            btestOut(i).F   = NaN;                    

        end

end

%% FOR GROUPED DATA
for i = 1 : length(X_grouped)

    if handles.displayMidResults == 1
        disp(' '); disp(['    ', X_grouped{i}.header, ' (n = ', num2str(length(X_grouped{i}.paired(:,1))), ')']); disp('    ======= '); disp(' ')    
    end

    if handles.testStatX == 0
        if handles.displayMidResults == 1
            X_grouped{i}.paired       
            X_grouped{i}.nonPaired
        end
    end

    X_grouped{i}.dataSample

        %% do Bartlett's TEST
        if sw.W(1) < 1 && sw.W(2) < 1 && sw.W(3) < 1

            cd(handles.path.statFunc)                

            % first with Trujillo-Ortiz's subfunction
            alpha  = handles.bartlett_pThreshold;
            X      = X_grouped{i}.dataSample;               

            btestOutGrouped{i} = Btest(X,alpha);

            if handles.displayMidResults == 1
                disp(['     ', 'Bartlett K-squared test: [btest.m]'])       
                % disp(['         variance      :', num2str(btestOutGrouped{i}.var(1)), ' ', num2str(btestOutGrouped{i}.var(2))])
                disp(['         deg. freedom  :', num2str(btestOutGrouped{i}.v)])
                disp(['         p-probability :', num2str(btestOutGrouped{i}.P)])
                disp(['         chi^2         :', num2str(btestOutGrouped{i}.X2)])
                disp(['         F-statistic   :', num2str(btestOutGrouped{i}.F)]);
                disp(' ');
            end
            cd(currDir)

        else

            btestOutGrouped{i}.v  = NaN;
            btestOutGrouped{i}.P  = NaN;
            btestOutGrouped{i}.X2 = NaN;
            btestOutGrouped{i}.F  = NaN;                

        end          


        %% If Data are normal (Gaussian) and deviations are homogeneous then you can do an ANOVA 
        % --> significant if p < 0.05

        if sw.W(1) < 1 && sw.W(2) < 1 && sw.W(2) < 3 && btestOutGrouped{i}.P > handles.bartlett_pThreshold                

            % The ANOVA test makes the following assumptions about the data in X:
                % All sample populations are normally distributed.
                % All sample populations have equal variance.
                % All observations are mutually independent.                



            % the Matlab Statistics Toolbox                
            alpha  = handles.anova_pThreshold;                
            X      = X_grouped{i}.paired;
            [anovaOut{i}.p, anovaOut{i}.table, anovaOut{i}.stats] = anova1(X, [], 'off');

                % Notches in the boxplot provide a test of group medians 
                % (see boxplot) different from the F test for means in
                % the ANOVA table             '
                if handles.displayMidResults == 1
                    disp(['     ', '1-way ANOVA: [anova1]'])                    
                    disp(['         p :', num2str(anovaOut{i}.p)])
                    disp(' ');
                end

        else
            if handles.displayMidResults == 1
                disp(['     WARNING: Conditions for ANOVA are not met!']);  disp(' ');
            end
            anovaOut{i}.p     = NaN;
            anovaOut{i}.table = NaN;
            anovaOut{i}.stats = NaN;

        end


        %% 4) Kruskal-Wallis test
            X      = X_grouped{i}.paired;
            [kruskWillis{i}.p, kruskWillis{i}.table, kruskWillis{i}.stats] = kruskalwallis(X, [], 'off');

            if handles.displayMidResults == 1
                disp(['     ', 'Kruskal-Wallis test: [kruskalwallis]'])                    
                disp(['         p :', num2str(kruskWillis{i}.p)])
                disp(' ');
            end


        %% If not you should go for a non-parametric test like : Wilcoxon-Mann-Whithney Test, unpaired sample
        % --> significant if p < 0.05

            % We use Giuseppe Cardillo's code for this, same as RANKSUM
            % in MATLAB,  % http://www.mathworks.com/matlabcentral/fileexchange/25830
            cd(handles.path.statFunc)
            doPaired = handles.wilcox_doPaired;

            if i ~=  length(X_grouped) % when only two columns / groups are to be analyzed
                if doPaired == 1
                    x1 = X_grouped{i}.nonPaired(~isnan(X_grouped{i}.nonPaired(:,1)),1);
                    x2 = X_grouped{i}.nonPaired(~isnan(X_grouped{i}.nonPaired(:,3)),2);                
                else
                    x1 = X_grouped{i}.paired(~isnan(X_grouped{i}.paired(:,1)),1);
                    x2 = X_grouped{i}.paired(~isnan(X_grouped{i}.paired(:,2)),2);
                end

                dispFullResultsFlag = 0;
                mww_Stats{i} = mwwtest(x1,x2,dispFullResultsFlag);                                

                if handles.displayMidResults == 1
                    disp(['     ', 'Wilcoxon-Mann-Whitney Test: [mwwtest.m], paired = ', num2str(doPaired)])                    
                    disp(['         method :', mww_Stats{i}.method])
                    disp(['         T      :', num2str(mww_Stats{i}.T)])
                    disp(['         U      :', num2str(mww_Stats{i}.U)])
                    % disp(['         mean   :', num2str(mww_Stats{i}.mean)])
                    % disp(['         std    :', num2str(mww_Stats{i}.std)])
                    % disp(['         z      :', num2str(mww_Stats{i}.z)])
                    % disp(['         p      :', num2str(mww_Stats{i}.p)])
                    disp(' ');                
                end

            % In mwwtest the method is determined automatically, so we
            % use the same method for the ransum -function as well
            if strcmp(mww_Stats{i}.method, 'Normal approximation') == 1
                handles.wilcox_method = 'approximate';
            else
                handles.wilcox_method = 'exact';
            end

            [rankSum{i}.p, rankSum{i}.h, rankSum{i}.stats] = ranksum(x1,x2,'alpha',handles.wilcox_pThreshold,'method',handles.wilcox_method);

            if handles.displayMidResults == 1
                disp(['     ', 'Wilcoxon-Mann-Whitney Test: [ranksum], paired = ', num2str(doPaired)'])                    
                disp(['         p       :', num2str(rankSum{i}.p)])
                disp(['         h       :', num2str(rankSum{i}.h)])
                % disp(['         zval    :', num2str(rankSum{i}.stats.zval)])
                % disp(['         ranksum :', num2str(rankSum{i}.stats.ranksum)])
                disp(' ');
            end

            else % for the situation for all the groups, pairwise matching not really possible

                mww_Stats{i} = [];                
                rankSum{i}   = [];

            end           

        %% Student's T test
        if i ~=  length(X_grouped) % when only two columns / groups are to be analyzed
            alpha  = handles.student_pThreshold; 
            x = X_grouped{i}.nonPaired(~isnan(X_grouped{i}.nonPaired(:,1)),1);
            y = X_grouped{i}.nonPaired(~isnan(X_grouped{i}.nonPaired(:,2)),2);
            [stT{i}.h, stT{i}.p, stT{i}.ci, stT{i}.stats] = ttest2(x,y,alpha);

            if handles.displayMidResults == 1
                disp(['     ', 'Student Two-sample t-test: [ttest2]'])                    
                disp(['         h       :', num2str(stT{i}.h)])
                disp(['         p       :', num2str(stT{i}.p)])                    
                % .ci      confidence interval
                % .stats   stats structure that can be passed easily to
                %          other statistical functions
                disp(' ');
            end

        else
            stT{i}.h = [];
            stT{i}.p = [];

        end
                
            

   end

   cd(currDir)
   %% OUTPUT STRUCTURE
   statsAdv.groupwise.shapWilk  = sw;
   statsAdv.groupwise.bartlett  = btestOut;
   statsAdv.bartlett            = btestOutGrouped;
   statsAdv.anova1              = anovaOut;
   statsAdv.mww_Stats           = mww_Stats;
   statsAdv.rankSum             = rankSum;
   statsAdv.kruskWillis         = kruskWillis;
   statsAdv.studentT            = stT;