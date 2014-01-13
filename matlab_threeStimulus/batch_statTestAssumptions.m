function stat_assumptions = batch_statTestAssumptions(inputVector, inputMatrix, conditions, condition, session, parameters, handles)

    %% DEBUG
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction    
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempStatAssumptions.mat';
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
    
    n = length(inputVector);
    parameters.stats.displayMidResults = 1;
    

    %% 1) You should test for the normality of distribution of the data. (Shapiro-Wilk normality test) 
    % ---> ok if p > 0.05

        % IMPLEMENTATION by: Ahmed Ben Saïda
        % http://www.mathworks.com/matlabcentral/fileexchange/13964
        % Shapiro-Wilk parametric hypothesis test of composite normality, for sample size 3<= n <= 5000. 
        % Based on Royston R94 algorithm. 
        % This test also performs the Shapiro-Francia normality test for platykurtic samples.        

        % inputs
        X = inputVector;
        alpha = parameters.stats.shapWilk_pThreshold;
        tail = 1; % default value, check swtest.m for explations

        % call the subfunction
        if n >= 3
            try
                [sw.H, sw.pValue, sw.W] = swtest(X, alpha, tail);
            catch err
                % err
                sw.H = NaN;
                sw.pValue = NaN;
                sw.W = NaN;                
            end
        else
            sw.H = NaN;
            sw.pValue = NaN;
            sw.W = NaN;
        end

        if parameters.stats.displayMidResults == 1
            disp(['     ', 'Shapiro-Wilk: [swtest.m]'])
            disp(['         H :', num2str(sw.H)])
            disp(['         p :', num2str(sw.pValue)])
            disp(['         W :', num2str(sw.W)])
        end

            % NULL hypothesis is that the distribution is NOT Gaussian so
            % if the p-value is higher than your alpha, then the
            % distribution is Gaussian, i.e. if W < 1
            if sw.W < 1
                stat_assumptions.shapiroWilk.isGaussian = 1;
            else
                stat_assumptions.shapiroWilk.isGaussian = 0;
            end
            stat_assumptions.shapiroWilk.W = sw.W;
            stat_assumptions.shapiroWilk.H = sw.H;
            stat_assumptions.shapiroWilk.p = sw.pValue;
        

    %% 2) If valid you test homogeneity of distribution of the deviation. (Bartlett's K-squared test) 
    % --> ok if p > 0.05

        if sw.W < 1 % if the distribution is indeed GAUSSIAN           

            % Use the subfunction provided by Antonio
            % Trujillo-Ortiz to reduce dependency on additional
            % toolboxes such as Statistics Toolbox that has the
            % function "barttest" (?)
            % http://www.mathworks.com/matlabcentral/fileexchange/3314-btest                                       
            data   = inputVector;
            sample = ones(length(data),1);
            X      = [data sample];
            alpha  = parameters.stats.bartlett_pThreshold;

            % first with Trujillo-Ortiz's subfunction
            stat_assumptions.bartlett = Btest(X,alpha);
            
            
            % DO TOTAL AS WELL (combining all the conditions: dark / dim
            % / bright)
            data = squeeze(inputMatrix(:, session, :));        
            sampleVector = (1 : 1 : size(data,1))';
            sample = repmat(sampleVector, 1, size(data,2));
            
            % vectorize
            data = data(:);
            sample = sample(:);
            X      = [data sample];
                        
            stat_assumptions.bartlett_total = Btest(X,alpha);
            

        if parameters.stats.displayMidResults == 1
            disp(['     ', 'Bartlett K-squared test (subgroup vectors, condition = ', conditions{condition}, '): [btest.m]'])                    
            disp(['         variance      :', num2str(stat_assumptions.bartlett.var)])
            disp(['         deg. freedom  :', num2str(stat_assumptions.bartlett.v)])
            disp(['         p-probability :', num2str(stat_assumptions.bartlett.P)])
            disp(['         chi^2         :', num2str(stat_assumptions.bartlett.X2)])
            disp(['         F-statistic   :', num2str(stat_assumptions.bartlett.F)]); 
            disp(['     ', 'Bartlett K-squared test (total, all conditions): [btest.m]'])                       
            fprintf(['         variance      :'])
            for i = 1 : length(conditions)
                fprintf([' ', num2str(stat_assumptions.bartlett_total.var(i))])
            end
            fprintf('\n')
            disp(['         deg. freedom  :', num2str(stat_assumptions.bartlett_total.v)])
            disp(['         p-probability :', num2str(stat_assumptions.bartlett_total.P)])
            disp(['         chi^2         :', num2str(stat_assumptions.bartlett_total.X2)])
            disp(['         F-statistic   :', num2str(stat_assumptions.bartlett_total.F)]); 
            disp(' ');
        end
            
            stat_assumptions.bartlett.var = NaN;
            stat_assumptions.bartlett.v   = NaN;
            stat_assumptions.bartlett.P   = NaN;
            stat_assumptions.bartlett.X2  = NaN;
            stat_assumptions.bartlett.F   = NaN;                    

        end

    %% 3) Sphericity-test: Mauchly's W
    
        % We can use the implementation "Huynh-Feldt epsilon general
        % procedure" by Matt Nelson to get the results of Mauchly's W test
        % http://www.mathworks.com/matlabcentral/fileexchange/22870-huynh-feldt-epsilon-general-procedure    
        %[EpsHF, EpsList, EpsGG, Mau] = GenCalcHFEps(inputVector, [], WInFacs)
        
        % Or the Mauchy implementation by Antonio Trujillo-Ortiz
        % NOTE! That some users on comments had noted with discrepancies
        % with the results obtained with this and the one from SPSS
        % http://www.mathworks.com/matlabcentral/fileexchange/3694-sphertest
        try
            stat_assumptions.mauchly.isTenable = Mauspher(inputVector,alpha) 
        catch
            stat_assumptions.mauchly.isTenable = NaN;
        end

         %{
         if parameters.stats.displayMidResults == 1
            disp(['     ', 'Mauchly''s W test ([GenCalcHFEps.m]'])
            disp(['         w      :', num2str(Mau.w)])
            disp(['         deg. freedom  :', num2str(Mau.df)])
            disp(['         p-probability :', num2str(Mau.p)])
            disp(['         chi^2         :', num2str(Mau.chiSq)])
         end
        %}

