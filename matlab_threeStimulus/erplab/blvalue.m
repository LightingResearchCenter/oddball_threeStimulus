%---------------------------------------------------------------------------------------------------
%-----------------base line mean value--------------------------------------------------------------
function blv = blvalue(ERP, chan, bin, blcorr)

%
% Baseline assessment
%
if ischar(blcorr)
      
      if ~strcmpi(blcorr,'no') && ~strcmpi(blcorr,'none')
            
            if strcmpi(blcorr,'pre')
                  bb = find(ERP.times==0);    % zero-time locked
                  aa = 1;
            elseif strcmpi(blcorr,'post')
                  bb = length(ERP.times);
                  aa = find(ERP.times==0);
            elseif strcmpi(blcorr,'all') || strcmpi(blcorr,'whole')
                  bb = length(ERP.times);
                  aa = 1;
            else
                  blcnum = str2num(blcorr); % in ms
                  
                  %
                  % Check & fix baseline range
                  %
                  if blcnum(1)<ERP.xmin*1000
                        blcnum(1) = ERP.xmin*1000; %ms
                  end
                  if blcnum(2)>ERP.xmax*1000
                        blcnum(2) = ERP.xmax*1000; %ms
                  end
                  
                  [xxx, cindex] = closest(ERP.times, blcnum); % 04/21/2011
                  aa = cindex(1); % ms to sample pos
                  bb = cindex(2); % ms to sample pos
            end
            blv = mean(ERP.bindata(chan,aa:bb, bin));
      else
            blv = 0;
      end
else      
      %
      % Check & fix baseline range
      %
      if blcorr(1)<ERP.xmin*1000
            blcorr(1) = ERP.xmin*1000; %ms
      end
      if blcorr(2)>ERP.xmax*1000
            blcorr(2) = ERP.xmax*1000; %ms
      end
      [xxx, cindex] = closest(ERP.times, blcorr); % 04/21/2011
      aa = cindex(1); % ms to sample pos
      bb = cindex(2); % ms to sample pos
      blv = mean(ERP.bindata(chan,aa:bb, bin));
end
