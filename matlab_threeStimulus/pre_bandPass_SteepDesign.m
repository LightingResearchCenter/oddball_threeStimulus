function [dataOut2] = pre_bandPass_SteepDesign(dataOut, lowPass, highPass, Fs, N, parameters, handles)

    debugMatFileName = 'tempPreBandbassFilterSteep.mat';
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

    % http://www.mathworks.com/help/dsp/ref/fdesign.bandpass.html

    % Design a constrained-band FIR equiripple filter 
    % Both stopband attenuation values are constrained to 60 dB
    
    fStopBand_lo = 0.1 * lowPass;
    fStopBand_hi = 4 * highPass;
    N = 36;
    
    d = fdesign.bandpass('N,Fst1,Fp1,Fp2,Fst2,C',N,fStopBand_lo,lowPass,highPass,fStopBand_hi,Fs);
    d.Stopband1Constrained = true; d.Astop1 = 1;
    d.Stopband2Constrained = true; d.Astop2 = 60;
    Hd = design(d,'equiripple');
    fvtool(Hd);
    measure(Hd)
    
    
    % continue at some point if needed, not functional really