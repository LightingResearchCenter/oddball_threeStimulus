function [ITPC, ITLC, ITLCN, ERSP, avWT, WTav, Induced] = analyze_wavelet_derivedParameters(ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi, noOfEpochs, parameters, handles)

    % Inter-Trial Phase Coherence (ITPC)
    ITPC = ITPC / noOfEpochs;

    % Inter-trial Linear Coherence (ITLC)
    ITLC = 1 / sqrt(noOfEpochs)*ITLC./sqrt(ITLCN);

    % Event-Related Spectral Perturbation (ERSP)
    ERSP = ERSP / noOfEpochs;

    % ?
    avWT = avWT / noOfEpochs;

    % ?
    WTav = WTav / noOfEpochs;

    % ?
    Induced = (WTavi-abs(avWTi)) / noOfEpochs;
        