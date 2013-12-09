
function ERP = checkerpzerolat(ERP)

if ~iserpstruct(ERP)
      fprintf('\nWARNING: checkeegzerolat() only works with ERP structure. This call was ignored.\n')
      return
end

auxtimes  = ERP.times;
[v indx]  = min(abs(auxtimes));
ERP.times = auxtimes - auxtimes(indx);
ERP.xmin  = min(ERP.times)/1000;
ERP.xmax  = max(ERP.times)/1000;
ERP.srate = round(ERP.srate);
if ERP.times(1)~=auxtimes(1)
      msg = ['\nWarning: zero time-locked stimulus latency values were not found.\n'...
      'Therefore, ERPLAB adjusted latency values at ERP.times, ERP.xmin,and ERP.xmax.\n\n'];
      fprintf(msg);
      fprintf('Time range is now [%.3f  %.3f] sec.\n', ERP.xmin, ERP.xmax )
else
      %fprintf('Zero latencies OK.\n')
end