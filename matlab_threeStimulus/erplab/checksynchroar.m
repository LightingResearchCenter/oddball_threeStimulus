% Author: Javier Lopez-Calderon & Steven Luck
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2012

function sync_status = checksynchroar(EEG)

sync_status = 1; % status for synchro. 1 means syncho; 0 means unsynchro

if ~isempty(EEG.reject.rejmanual)
      if length(EEG.reject.rejmanual) ~= max([EEG.EVENTLIST.eventinfo.bepoch])
            sync_status = 2; % deleted epochs detected
            return
      end
end

fprintf('\n---------------------------------------------------------\n');
fprintf('Testing artifact info synchronization I: EEG.reject.rejmanual vs EEG.epoch.eventflag...\n');
fprintf('---------------------------------------------------------\n\n');

nepoch = EEG.trials;
for i=1:nepoch
      cflag = EEG.epoch(i).eventflag; % flag(s) from event(s) within this epoch
      if iscell(cflag)
            %cflag = cell2mat(cflag);
            cflag = uint16([cflag{:}]); % giving some problems with uint16 type of flags
      end
      laten = EEG.epoch(i).eventlatency;% latency(ies) from event(s) within this epoch
      if iscell(laten)
            laten = cell2mat(laten);
      end
      
      indxtimelock = find(laten == 0,1,'first'); % catch zero-time locked code position,
      flag  = cflag(indxtimelock);
      
      if ~isempty(EEG.reject.rejmanual)
            if flag>0 && flag<=255  && EEG.reject.rejmanual(i)==0; %
                  sync_status = 0;
                  iflag = find(bitget(flag,1:8));
                  fprintf('Epoch # %g is not marked as artifactual but flag(s) # %s is(are) set.\n', i, num2str(iflag));
            elseif flag==0 && EEG.reject.rejmanual(i)==1
                  sync_status = 0;
                  fprintf('Epoch # %g is marked as artifactual but no flag is set.\n', i);
            end
      else
            if flag>0 && flag<=255 %
                  sync_status = 0;
                  iflag = find(bitget(flag,1:8));
                  fprintf('Epoch # %g is not marked as artifactual but flag(s) # %s is(are) set.\n', i, num2str(iflag));
            end
      end
end
if sync_status
      disp('Ok!')
else
      fprintf('\nFail!\n\n');
end

fprintf('\n---------------------------------------------------------\n');
fprintf('Testing artifact info synchronization II: EEG.reject.rejmanual vs EEG.EVENTLIST.eventinfo.flag...\n');
fprintf('---------------------------------------------------------\n\n');

nitem = length(EEG.EVENTLIST.eventinfo);
for i=1:nitem
      flag   = EEG.EVENTLIST.eventinfo(i).flag;
      bepoch = EEG.EVENTLIST.eventinfo(i).bepoch;
      if bepoch>0
            if ~isempty(EEG.reject.rejmanual)
                  if flag>0 && flag<=255 && EEG.reject.rejmanual(bepoch) == 0;
                        sync_status = 0;
                        iflag = find(bitget(flag,1:8));
                        fprintf('Epoch # %g is not marked as artifactual but flag(s) # %s is(are) set.\n', i, num2str(iflag));
                  elseif flag==0 && EEG.reject.rejmanual(bepoch)==1
                        sync_status = 0;
                        fprintf('Item # %g is not marked as artifactual but its corresponding epoch # %g.\n', i, bepoch);
                  end
            else
                  if flag>0 && flag<=255
                        sync_status = 0;
                        iflag = find(bitget(flag,1:8));
                        fprintf('Epoch # %g is not marked as artifactual but flag(s) # %s is(are) set.\n', i, num2str(iflag));
                  end
            end
      end
end
if sync_status
      disp('Ok!')
else
      fprintf('\nFail!\n\n');
end
return