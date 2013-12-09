
function EEG = checkeegzerolat(EEG)

if ~iseegstruct(EEG)
      fprintf('\nWARNING: checkeegzerolat() only works with EEG structure. This call was ignored.\n')
      return
end
if isempty(EEG.epoch)
      fprintf('\nWARNING: checkeegzerolat() only works for epoched datasets. This call was ignored.\n')
      return
end
nepoch = EEG.trials;
for i=1:nepoch
      latcell = EEG.epoch(i).eventlatency;
      if iscell(latcell)
            latcell  =  cell2mat(latcell);
            [v indx] = min(abs(latcell));
            EEG.epoch(i).eventlatency = num2cell(latcell - latcell(indx));
      else
            [v indx] = min(abs(latcell));
            EEG.epoch(i).eventlatency = latcell - latcell(indx);
      end
end
auxtimes  = EEG.times;
[v indx]  = min(abs(auxtimes));
EEG.times = auxtimes - auxtimes(indx);
EEG.xmin  = min(EEG.times)/1000;
EEG.xmax  = max(EEG.times)/1000;
EEG.srate = round(EEG.srate);

EEG = eeg_checkset( EEG );

if EEG.times(1)~=auxtimes(1)
      msg = ['\nWarning: zero time-locked stimulus latency values were not found.\n'...
      'Therefore, ERPLAB adjusted latency values at EEG.epoch.eventlatency, EEG.times, EEG.xmin,and EEG.xmax.\n\n'];
      fprintf(msg);
      fprintf('Time range is now [%.3f  %.3f] sec.\n', EEG.xmin, EEG.xmax )
else
      fprintf('Zero latencies OK.\n')
end