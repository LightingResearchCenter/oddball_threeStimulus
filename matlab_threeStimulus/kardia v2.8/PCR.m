%% ---- PCR
function [HRmean,HRbsl]=PCR(t,TS,T0,T1,step,unit)

    % get heart period
    [hp,t]=ecg_hp(t,'instantaneous');

    % calculate baseline HR
    HRbsl=ecg_stat(t,T0,TS,'mean',unit); % call ecg_stat

    % calculate mean HR changes in variable window sizes defined by step
    nboxes=floor(T1/step); % number of boxes that fit in analysis window
    HRmean=[];
    count=0;
    for i=1:nboxes
        mhr=ecg_stat(t,TS+count*step,TS+step*i,'mean',unit);
        count=count+1;
        HRmean=[HRmean mhr];
    end
end