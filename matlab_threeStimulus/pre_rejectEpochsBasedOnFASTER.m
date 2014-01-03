function epochsOut = pre_rejectEpochsBasedOnFASTER(epochsIn, artifactIndices, stimType, erpType, parameters, handles)
    
    debugAverageHere = 0;    

    %{
    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempRejectEpochs.mat';
        if nargin == 0
            load('debugPath.mat')
            load(fullfile(path.debugMATs, debugMatFileName))
            close all
            debugAverageHere = 1;
        else
            if handles.flags.saveDebugMATs == 1
                if 1 == 1 % ~strcmp(erpType, 'standard')
                    % do not save for standard tone as there are so many
                    % trials that debugging and developing of this function
                    % is so much slower compared to target and distracter
                    path = handles.path;
                    save('debugPath.mat', 'path')
                    save(fullfile(path.debugMATs, debugMatFileName))            
                end        
            end
        end 
    end    
    stimType
    erpType    
    artifactIndices
    %}

    epochsOut = epochsIn;
    noOfEpochsIn = length(epochsIn.ERP);
    
    for ep = 1 : noOfEpochsIn                
        
        artifactYes = artifactIndices(ep,:);                     
        % subplot(2,1,1); plot(epochsOut.ERP{ep}(:,1:4))
        
        epochsOut.ERP{ep}(:, artifactYes == 1) = NaN;        
        % subplot(2,1,2);plot(epochsOut.ERP{ep}(:,1:4));title(num2str(artifactYes)); pause(1.5)
                
        
    end
    
    %{
    if debugAverageHere == 1
       
        for ep = 1 : noOfEpochsIn
            epochAverage(ep,:,:) = epochsOut.ERP{ep}(:, 1:4);
        end
        
        meanDebug = squeeze(nanmean(epochAverage,1));
        plot(meanDebug); 
        notNans = sum(~artifactIndices)
        title(['n = ', num2str(notNans)])
        whos
        
    end
    %}
    
    