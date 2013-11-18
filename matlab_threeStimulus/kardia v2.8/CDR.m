%% ---- CDR
function [baselineHR,medians,Vsec]=CDR(sig,TS,baseline)

    % EXAMPLE USE from kardia.m
    % [HRbsl,HR_mean,values]=CDR(Rdata, lat, 15);
    % Rdata  -  
    whos
        
    % transpose data vector
    if size(sig,1)>size(sig,2)
        sig=sig';
    end

    % create heart period vector
    [hp,t]=ecg_hp(sig,'instantaneous');
    whos

    % define T0 to calculate baseline heart rate
    T0=TS-baseline;
    baselineHR=ecg_stat(t,T0,TS,'mean','bpm'); % baseline Heart Rate

    % calculate mean heart rate sec by sec
    Vsec=[];
    count=0;
    for i=1:80
        V=ecg_stat(t,TS+count,TS+i,'mean','bpm');
        count=count+1;
        Vsec=[Vsec V];
    end

    % calculate medians
    medians=[median(Vsec(1:3)) median(Vsec(4:6)) median(Vsec(7:11)) median(Vsec(12:16))...
        median(Vsec(17:23)) median(Vsec(24:30)) median(Vsec(31:37)) median(Vsec(38:50))...
        median(Vsec(51:63)) median(Vsec(64:76))];
end