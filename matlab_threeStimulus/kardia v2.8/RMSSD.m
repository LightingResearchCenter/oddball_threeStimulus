%% ---- RMSSD
function y=RMSSD(sig)

    dsig=diff(sig);
    y=dsig.^2;
    y=sqrt(mean(y));
end