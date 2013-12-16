% MAIN FUNCTION for Learning Oddball EEG ANALYSIS
function MAIN_ANALYSIS()



    % Petteri Teikari, petteri.teikari@gmail.com, 2013
    % Lighting Research Center, Rensselaer Polytechnic Institute, Troy, NY
    close all
    clear all

    % http://www.network-science.de/ascii/, font "rectangle"
    
        cl = fix(clock); hours = num2str(cl(4)); % get the current time
        if cl(5) < 10; mins = ['0', num2str(cl(5))]; else mins = num2str(cl(5)); end
        disp(' ');
        disp('       _   _ _       _ _    _____         _         _         ')
        disp(' ___ _| |_| | |_ ___| | |  |  _  |___ ___| |_ _ ___|_|___     ')
        disp('| . | . | . | . | .''| | |  |     |   | .''| | | |_ -| |_ -|  ')
        disp('|___|___|___|___|__,|_|_|  |__|__|_|_|__,|_|_  |___|_|___|    ')
        disp('                                           |___|              ')  
        disp(' ')
        disp(['Initiated: ', date, ', ', hours, ':', mins])
        disp('-------'); 
    
    %% General Settings
    
        % i.e. like where the folders are, fonts to be used, colors, etc.
        [handles, ~] = init_DefaultSettings(); % use a subfunction        

    %% Parameters for ANALYSIS
    
        % i.e. like artifact rejection threolds, filter cutoffs,
        % numerical parameters for EEG analysis
        handles.parameters = init_DefaultParameters(handles); % use a subfunction        
    
    %% Define input
    
        %{
        filesIn = {'je_01_40.bdf'; 'je_02_40.bdf'; 'je_03_40.bdf'; 'je_04_40.bdf'; 'je_01_0.bdf'; 'je_02_0.bdf'; 'je_03_0.bdf'; 'je_04_0.bdf'; 'je_01_10.bdf'; 'je_02_10.bdf'; 'je_03_10.bdf'; 'je_04_10.bdf'; ...
                   'bl_01_10.bdf'; 'bl_02_10.bdf'; ...
                   'bl_03_10.bdf'; 'bl_04_10.bdf'; 'bl_01_40.bdf'; 'bl_02_40.bdf'; 'bl_03_40.bdf'; 'bl_04_40.bdf'; 'bl_01_0.bdf'; 'bl_02_0.bdf'; 'bl_03_0.bdf'; 'bl_04_0.bdf'; ...
                   'ka_01_10.bdf'; 'ka_02_10.bdf'; 'ka_03_10.bdf'; 'ka_04_10.bdf'; 'ka_01_40.bdf'; 'ka_02_40.bdf'; 'ka_03_40.bdf'; 'ka_04_40.bdf'; 'ka_01_0.bdf'; 'ka_02_0.bdf'; 'ka_03_0.bdf'; 'ka_04_0.bdf'; ...
                   'da_01_10.bdf'; 'da_02_10.bdf'; 'da_03_10.bdf'; 'da_04_10.bdf'; 'da_01_40.bdf'; 'da_02_40.bdf'; 'da_03_40.bdf'; 'da_04_40.bdf'; 'da_01_0.bdf'; 'da_02_0.bdf'; 'da_03_0.bdf'; 'da_04_0.bdf'; ...
                   'sh_01_40.bdf'; 'sh_02_40.bdf'; 'sh_03_40.bdf'; 'sh_04_40.bdf'; 'sh_01_0.bdf'; 'sh_02_0.bdf'; 'sh_03_0.bdf'; 'sh_04_0.bdf'; 'sh_01_10.bdf'; 'sh_02_10.bdf'; 'sh_03_10.bdf'; 'sh_04_10.bdf'; ...
                   'dk_01_10.bdf'; 'dk_02_10.bdf'; 'dk_03_10.bdf'; 'dk_04_10.bdf'; 'dk_01_40.bdf'; 'dk_02_40.bdf'; 'dk_03_40.bdf'; 'dk_04_40.bdf'; 'dk_01_0.bdf'; 'dk_02_0.bdf'; 'dk_03_0.bdf'; 'dk_04_0.bdf'; ...        
                   'jo_01_40.bdf'; 'jo_02_40.bdf'; 'jo_03_40.bdf'; 'jo_04_40.bdf'; 'jo_01_0.bdf'; 'jo_02_0.bdf'; 'jo_03_0.bdf'; 'jo_04_0.bdf'; 'jo_01_10.bdf'; 'jo_02_10.bdf'; 'jo_03_10.bdf'; 'jo_04_10.bdf'; ...
                   'cr_01_10.bdf'; 'cr_02_10.bdf'; 'cr_03_10.bdf'; 'cr_04_10.bdf'; 'cr_01_40.bdf'; 'cr_02_40.bdf'; 'cr_03_40.bdf'; 'cr_04_40.bdf'; 'cr_01_0.bdf'; 'cr_02_0.bdf'; 'cr_03_0.bdf'; 'cr_04_0.bdf'; ...
                    'hh_01_40.bdf'; 'hh_02_40.bdf'; 'hh_03_40.bdf'; 'hh_04_40.bdf'; 'hh_01_0.bdf'; 'hh_02_0.bdf'; 'hh_03_0.bdf'; 'hh_04_0.bdf'; 'hh_01_10.bdf'; 'hh_02_10.bdf'; 'hh_03_10.bdf'; 'hh_04_10.bdf'; ...
                   'ha_01_10.bdf'; 'ha_02_10.bdf'; 'ha_03_10.bdf'; 'ha_04_10.bdf'; 'ha_01_40.bdf'; 'ha_02_40.bdf'; 'ha_03_40.bdf'; 'ha_04_40.bdf'; 'ha_01_0.bdf'; 'ha_02_0.bdf'; 'ha_03_0.bdf'; 'ha_04_0.bdf'};
        %}     
     
       filesIn = {'cr_03_0.bdf'; 'cr_04_0.bdf'; ...
                    'hh_01_40.bdf'; 'hh_02_40.bdf'; 'hh_03_40.bdf'; 'hh_04_40.bdf'; 'hh_01_0.bdf'; 'hh_02_0.bdf'; 'hh_03_0.bdf'; 'hh_04_0.bdf'; 'hh_01_10.bdf'; 'hh_02_10.bdf'; 'hh_03_10.bdf'; 'hh_04_10.bdf'; ...
                   'ha_01_10.bdf'; 'ha_02_10.bdf'; 'ha_03_10.bdf'; 'ha_04_10.bdf'; 'ha_01_40.bdf'; 'ha_02_40.bdf'; 'ha_03_40.bdf'; 'ha_04_40.bdf'; 'ha_01_0.bdf'; 'ha_02_0.bdf'; 'ha_03_0.bdf'; 'ha_04_0.bdf'};

        % filesIn = {'bl_02_40.bdf'};     
        fileNameIn = filesIn;
        
        handles.inputFile = fileNameIn;
        inputFiles = fullfile(handles.path.dataFolder, fileNameIn);
        if ~iscell(inputFiles) % for only one file
            numberOfFiles = 1;
        else
            numberOfFiles = length(inputFiles);
        end
           
        
    %% GO THROUGH THE FILES    
    
        hiFreq = [2 6 10 14 16 20 24 28 32 36 40];
        loFreq = [0.001 0.01 0.1 1];        
        scales = 3:1:10;
        % cnvFreq = [1 2 3 4 5 6 7 8 9 10 12 14 16];
        
        hiFreq = [14 20];
        loFreq = [0.01];        
        scales = [8 9];
        
        % only one value
        
        hiFreq = [handles.parameters.filter.bandPass_ERP_hiFreq];
        loFreq = [handles.parameters.filter.bandPass_ERP_loFreq];
        scales = [handles.parameters.ep_den.scales_postStim];
                
        numberOfCombinations = length(hiFreq) * length(loFreq) * length(scales)          
                
        for i = 1 : length(filesIn)
            
            tic
            close all            
            disp(['file: ', num2str(i), ' / ', num2str(length(filesIn))])
        
            for ij = 1 : length(loFreq)
                
                for ik = 1 : length(hiFreq)
                 
                    for il = 1 : length(scales)
                    
                        handles.inputFile = filesIn{i}; 
                        inputFiles = fullfile(handles.path.dataFolder, filesIn{i});
                        % disp(['cut-off high freq = ', num2str(hiFreq(i))])
                        handles.parameters.ep_den.scales_postStim = scales(il);
                        handles.parameters.filter.bandPass_ERP_hiFreq = hiFreq(ik);
                        handles.parameters.filter.bandPass_ERP_loFreq = loFreq(ij);

                        %% IMPORT THE DATA
                        % You could modify what is inside here if you have EEG recorded
                        % with some other system than with BioSemi ActiveTwo
                        [dataMatrix, triggers, info, handles.parameters.EEG.srate] = IMPORT_eegData(i, fileNameIn, inputFiles, handles);

                        %% PROCESS the ERP
                        [epochs, analyzed, TF, dataMatrix_filtGeneral, alpha, powers, handles] = ...
                            PROCESS_singleFile(inputFiles, fileNameIn, dataMatrix, triggers, handles);        
                            % Check later what is coming out in handles and try to refer to
                            % exact variables changed (like sampleRate)                      

                        %% PLOT
                        %{
                        PLOT_singleFile(epochs_target, epochs_distracter, epochs_standard, ...
                         epochsEP_target, epochsEP_distracter, epochsEP_standard, ...
                         analyzed_target, analyzed_distracter, analyzed_standard, ...
                         analyzed_aux, analyzed_extraSensors, analyzed_fractal, analyzed_TF, ...
                         dataMatrix_filtGeneral, alpha, powers, info, handles)
                         %}
                        
                    end
                end
            end

            tElapsed(i) = toc;
            disp(['TIME FOR COMPUTATION: ', num2str(tElapsed(i)), ' sec'])
            whos

        end           
        
        timePerFile = mean(tElapsed)
        timeTotal = sum(tElapsed)