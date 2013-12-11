function epochsOut = pre_rejectEpochsBasedOnFASTER(epochsIn, artifactIndices, parameters, handles)
    
    epochsOut = epochsIn;
    noOfEpochsIn = length(epochsIn.ERP);
    
    for ep = 1 : noOfEpochsIn                
        artifactYes = artifactIndices(ep,:);
        %subplot(2,1,1); plot(epochsOut.ERP{ep})
        epochsOut.ERP{ep}(:, artifactYes) = NaN;
        %subplot(2,1,2);plot(epochsOut.ERP{ep});title(num2str(artifactYes)); pause(1.5)
    end