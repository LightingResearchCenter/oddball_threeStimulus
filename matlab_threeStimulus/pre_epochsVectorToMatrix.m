function EEG_mat = pre_epochsVectorToMatrix(EEG,epochLength)

    [noOfChannels, noOfSamples] = size(EEG);
    noOfEpochs = noOfSamples / epochLength;

    EEG_mat = zeros(noOfChannels, epochLength, noOfEpochs);

    for ep = 1 : noOfEpochs

        ind1 = (ep-1)*epochLength + 1;
        ind2 = epochLength*ep;

        % disp([ind1 ind2])
        EEG_mat(:,:,ep) = EEG(:,ind1:ind2);

    end