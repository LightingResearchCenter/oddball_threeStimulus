function [ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi] = analyze_wavelet_derivedParametersAccum(WT, ITPC, ITLC, ITLCN, ERSP, avWT, WTav, avWTi, WTavi, parameters, handles)

    % Inter-Trial Phase Coherence (ITPC)
    ITPC = ITPC + WT./abs(WT);

    % Inter-trial Linear Coherence (ITLC)
    ITLC = ITLC+WT;    
    ITLCN = ITLCN+abs(WT).^2; % Inter-trial Linear Coherence (ITLCN) NORMALIZED?

    % Event-Related Spectral Perturbation (ERSP)
    ERSP = ERSP+abs(WT).^2;

    avWT = avWT+WT;

    WTav = WTav+abs(WT);

    % INDUCED
    avWTi = avWTi+WT;
    WTavi = WTavi+abs(WT);