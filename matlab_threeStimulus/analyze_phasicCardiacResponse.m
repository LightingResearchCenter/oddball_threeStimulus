function PCROut = analyze_phasicCardiacResponse(epochs_ECG, parameters, handles)

    % Use KARDIA
    % http://sourceforge.net/projects/mykardia/
    % http://dx.doi.org/10.1016/j.cmpb.2009.10.002
    PCROut = [];