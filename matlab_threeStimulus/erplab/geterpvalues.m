% Usage:
% >> VALUES  = geterpvalues(ERP, fname, latency, chanArray, op, dig)
%
% or
%
% >> [VALUES L] = geterpvalues(ERP, fname, latency, chanArray, op, dig)
%
% INPUTS
%
% ERP        - ERP structure
% fname      - name of text file for output. e.g. 'S03_N100_peak.txt'
% latency    - one or two latencies in msec. e.g. [80 120]
% chanArray  - index(es) of channel(s) you want to extract the information. e.g. [10 22  38 39 40]
% op         - option. Any of these:
%                'instabl'    finds the relative-to-baseline instantaneous value at the specified latency.
%                'peakampbl'  finds the relative-to-baseline peak value
%                             between two latencies. See polpeak and sampeak.
%                'peaklatbl'  finds peak latency between two latencies. See polpeak and sampeak.
%                'meanbl'     calculates the relative-to-baseline mean amplitude value between two latencies.
%                'area'       calculates the area under the curve value between two latencies.
%                '50arealat' calculates the latency corresponding to the 50% area sample between two latencies.
%                'areaz'      calculates the area under the curve value. Lower and upper limit of integration
%                             are automatically found starting with a seed
%                             latency.
% coi        - component of interest (1 or 2) (only for 'areaz')
% dig        - number of digit to use, for precision, used to write the text file for output. Default is 4
% polpeak    - polarity of peak:  1=maximum (default), 0=minimum
% sampeak    - number of points in the peak's neighborhood (one-side) (0 default)
%
%
% OUTPUTS
%
% VALUES     - matrix of values. bin(s) x channel(s). geterpvalues() use always all bins.
% L          - Latencies structure: fields are:
%              "value"  : latency in msec
%              "ilimit" : limits of integration in msec in case of using "area" or "areaz" as an option
%
% text file  - text file containing formated values.
%
% Author: Javier Lopez-Calderon & Steven Luck
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009

%b8d3721ed219e65100184c6b95db209bb8d3721ed219e65100184c6b95db209b
%
% ERPLAB Toolbox
% Copyright © 2007 The Regents of the University of California
% Created by Javier Lopez-Calderon and Steven Luck
% Center for Mind and Brain, University of California, Davis,
% javlopez@ucdavis.edu, sjluck@ucdavis.edu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout  = geterpvalues(ERP, latency, binArray, chanArray, op1, blc, coi, polpeak, sampeak, localopt, frac, fracmearep)

if nargin<1
      help geterpvalues
      return
end
if nargin<2
      error('ERROR geterpvalues(): You must specify ERP struct and latency(ies), at least.')
end
if nargin<12
      fracmearep = 0; %0=write a NaN when frac measure is not found.  1 = export frac absolute peak when frac local peak is not found.; 2=shows error message
end
if nargin<11
      frac = 0.5; % 50% area latency (if needed, by default)
end
if nargin<10
      localopt = 0; % 0=write a NaN when local peak is not found.  1=export absolute peak when local peak is not found.
end
if nargin<9
      sampeak = 0; % absolute peak. No neighbor samples
end
if nargin<8
      polpeak = 1; % positive
end
if nargin<7
      coi = 0; % 0= as it is; 1=first component; 2=2nd component
end
if nargin<6
      blc = 'pre';
end
if nargin<5
      op1 = 'instabl';
end
if nargin<4
      chanArray = 1:ERP.nchan;
end
if nargin<3
      binArray = 1:ERP.nbin;
end
if ischar(blc)
      
      blcnum = str2num(blc);
      
      if isempty(blcnum)
            if ~ismember(blc,{'no','none','pre','post','all','whole'})
                  msgboxText =  'Invalid baseline range dude!';
                  title        =  'ERPLAB: geterpvalues() baseline input';
                  errorfound(msgboxText, title);
                  return
            end
      else
            if size(blcnum,1)>1 || size(blcnum,2)>2
                  msgboxText =  'Invalid baseline range, dude!';
                  title        =  'ERPLAB: geterpvalues() baseline input';
                  errorfound(msgboxText, title);
                  return
            end
      end
end
if ~ismember({op1}, {'instabl', 'meanbl', 'peakampbl', 'peaklatbl', 'fpeaklat',...
            'areat', 'areap', 'arean','areazt','areazp','areazn','fareatlat',...
            'fareaplat', 'fninteglat', 'fareanlat', 'ninteg','nintegz' });
      msgboxText =  [op1 ' is not a valid option for geterpvalues!'];
      title = 'ERPLAB: geterpvalues wrong inputs';
      errorfound(msgboxText, title);
      return
end
if isempty(coi)
      coi = 1;
end
if nargout==1
      condf = 0; % only includes area values
elseif nargout==2
      condf = 1; % include latency values and limits...
else
      error('ERPLAB says: error at geterpvalues(). Too many output arguments!')
end

fs      = ERP.srate;
pnts    = ERP.pnts;
nbin    = length(binArray);
nchan   = length(chanArray);
nlat    = length(latency);
VALUES  = zeros(nbin,nchan);
LATENCY = struct([]);

% latsamp = [find(ERP.times>=latency(1), 1, 'first') find(ERP.times<=550, 1, 'last')];
% % its fields are "value" and "ilimits"
% toffsa  = round(ERP.xmin*fs);                    % in samples
% latsamp = round(latency*fs/1000) - toffsa + 1;   % msec to samples

[worklate{1:nbin,1:nchan}] = deal(latency); % specified latency(ies) for getting measurements.

if length(latency)==2
      
      %latsamp = [find(ERP.times>=latency(1), 1, 'first') find(ERP.times<=latency(2), 1, 'last')];
      [xxx, latsamp] = closest(ERP.times, latency);
      
      if latsamp(1)<-1
            msgboxText =  'ERROR: The onset of your latency window cannot be more than 2 samples away from the real onset ';
            tittle = 'ERPLAB: geterpvalues()';
            errorfound(msgboxText, tittle);
            varargout = {[]};
            return
      end
      if latsamp(2)>pnts+2
            msgboxText =  'ERROR: The offset of your latency window cannot be more than 2 samples away from the real offset';
            tittle = 'ERPLAB: geterpvalues()';
            errorfound(msgboxText, tittle);
            varargout = {[]};
            return
      end
      if latsamp(1)<1
            latsamp(1) = 1;
            fprintf('\n%s\n', repmat('*',1,60));
            fprintf('WARNING: Lower latency limit %.1f ms was adjusted to %.1f ms \n', latency(1), 1000*ERP.xmin);
            fprintf('%s\n\n', repmat('*',1,60));
      end
      if latsamp(2)>pnts
            latsamp(2) = pnts;
            fprintf('\n%s\n', repmat('*',1,60));
            fprintf('WARNING: Upper latency limit %.1f ms was adjusted to %.1f ms \n', latency(2), 1000*ERP.xmax);
            fprintf('%s\n\n', repmat('*',1,60));
      end
elseif length(latency)==1
      
      %latsamp = find(ERP.times>=latency(1), 1, 'first');
      [xxx, latsamp] = closest(ERP.times, latency(1));
      
      if latsamp(1)<-1 || latsamp(1)>pnts+2
            msgboxText{1} =  'ERROR: The specified latency is more than 2 samples away from the ERP window.';
            tittle = 'ERPLAB: geterpvalues()';
            errorfound(msgboxText, tittle);
            varargout = {[]};
            return
      end
      if latsamp(1)<1
            latsamp(1) = 1;
            fprintf('\n%s\n', repmat('*',1,60));
            fprintf('WARNING: Latency %.1f ms was adjusted to %.1f ms \n', latency(2), 1000*ERP.xmin);
            fprintf('%s\n\n', repmat('*',1,60));
      elseif latsamp(1)>pnts
            latsamp(1) = pnts;
            fprintf('\n%s\n', repmat('*',1,60));
            fprintf('WARNING: Latency %.1f ms was adjusted to %.1f ms \n', latency(2), 1000*ERP.xmax);
            fprintf('%s\n\n', repmat('*',1,60));
      end
else
      error('Wrong number of latencies...')
end

try      
      for b=1:nbin
            for ch = 1:nchan
                  if nlat==1   % 1 latency was specified
                        if strcmpi(op1,'areazt') || strcmpi(op1,'areazp') || strcmpi(op1,'areazn')
                              
                              %
                              % get area (automatic limits, 1 seed latency)
                              %
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;
                              
                              switch op1
                                    case 'areazt'
                                          op2 = 'autot';
                                    case 'areazp'
                                          op2 = 'autop';
                                    case 'areazn'
                                          op2 = 'auton';
                              end
                              
                              [A L il] =  areaerp(dataux, ERP.srate, latsamp, op2, coi);
                              worklate{b,ch} = ((il-1)/fs + ERP.xmin)*1000;   % integratin limits
                              VALUES(b,ch)  = A;
                              
                        elseif strcmpi(op1,'nintegz')
                              
                              %
                              % get numerical integration (automatic limits)
                              %
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;
                              
                              %if condf
                              [A L il]  =  areaerp(dataux, ERP.srate,latsamp, 'auto', coi) ;
                              worklate{b,ch} = ((il-1)/fs + ERP.xmin)*1000; % integratin limits
                              VALUES(b,ch)  = A;
                              
                        elseif strcmpi(op1,'instabl')
                              
                              %
                              % get instantaneous amplitud (at 1 latency)
                              %
                              blv  = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              VALUES(b,ch) = ERP.bindata(chanArray(ch),latsamp, binArray(b)) - blv;
                              
                        else
                              stre1  = 'ERROR in geterpvalues.m: You must enter 2 latencies for ';
                              stre2  = '''meanbl'', ''peakampbl'', ''peaklatbl'', or ''area''';
                              strerr = [stre1 stre2];
                              error( strerr )
                        end
                  else % between 2 latencies measurements.
                        
                        if strcmpi(op1,'meanbl')
                              
                              %
                              % get mean value
                              %
                              blv  = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value                              
                              VALUES(b,ch)  = mean(ERP.bindata(chanArray(ch),latsamp(1):latsamp(2), binArray(b))) - blv;
                              
                        elseif strcmpi(op1,'peakampbl')
                              
                              %
                              % get peak amplitude
                              %
                              blv  = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              
                              try
                                    dataux  = ERP.bindata(chanArray(ch),latsamp(1)-sampeak:latsamp(2)+sampeak, binArray(b)) - blv;  % no filtered
                                    timex   = ERP.times(latsamp(1)-sampeak:latsamp(2)+sampeak);
                                    %[valmaxf posmaxf ]= localpeak(dataux, sampeak, polpeak);
                              catch
                                    error(sprintf('ERPLAB says: The requested measurement range (%g - %g +/- %g ms) exceeds the time range of the data (%.1f - %.1f ms).',...
                                          round(latency(1)), round(latency(2)), round(1000*sampeak/fs), ERP.xmin*1000, ERP.xmax*1000))
                              end
                              
                              if localopt==1 %0=writes a NaN when local peak is not found.  1=export absolute peak when local peak is not found.
                                    localoptstr = 'abs';
                              else
                                    localoptstr = 'NaN';
                              end
                              
                              %%% if strcmpi(p.Results.Measure,'amplitude')
                              %%%       measure = vlocalpf;
                              %%% elseif strcmpi(p.Results.Measure,'peaklat')
                              %%%       measure = latlocalpeak;
                              %%% elseif strcmpi(p.Results.Measure,'fraclat')
                              %%%       measure = latfracpeak;
                              %%% end
                              
                              %%% [measure, latlocalpeak, latfracpeak, errorcode]
                              
                              [valx latpeak] = localpeak(dataux, timex, 'Neighborhood',sampeak, 'Peakpolarity', polpeak, 'Sfactor', 4, 'Measure','amplitude',...
                                    'Peakreplace', localoptstr);
                              
                              if isempty(valx)
                                    error('Peak-related measurement failed...')
                              end
                              
                              worklate{b,ch} = latpeak; %((il-1)/fs + ERP.xmin)*1000;
                              VALUES(b,ch)  = valx; % value of amplitude
                              
                        elseif strcmpi(op1,'peaklatbl')
                              
                              %
                              % get peak latency
                              %
                              blv  = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              
                              try
                                    dataux  = ERP.bindata(chanArray(ch),latsamp(1)-sampeak:latsamp(2)+sampeak, binArray(b)) - blv;  % no filtered
                                    timex   = ERP.times(latsamp(1)-sampeak:latsamp(2)+sampeak);
                                    %[valmaxf posmaxf ]= localpeak(dataux, sampeak, polpeak);
                              catch
                                    error(sprintf('ERPLAB says: The requested measurement range (%g - %g +/- %g ms) exceeds the time range of the data (%.1f - %.1f ms).',...
                                          round(latency(1)), round(latency(2)), round(1000*sampeak/fs), ERP.xmin*1000, ERP.xmax*1000))
                              end
                              
                              if localopt==1 %0=writes a NaN when local peak is not found.  1=export absolute peak when local peak is not found.
                                    localoptstr = 'abs';
                              else
                                    localoptstr = 'NaN';
                              end
                              
                              valx = localpeak(dataux, timex, 'Neighborhood', sampeak, 'Peakpolarity', polpeak, 'Sfactor', 4, 'Measure','peaklat',...
                                    'Peakreplace', localoptstr);
                              
                              if isempty(valx)
                                    error('Peak-related measurement failed...')
                              end
                              
                              VALUES(b,ch) = valx;
                              
                        elseif strcmpi(op1,'areat') || strcmpi(op1,'areap') || strcmpi(op1,'arean')
                              
                              %
                              % get area
                              %
                              switch op1
                                    case 'areat'
                                          op2 = 'total';
                                    case 'areap'
                                          op2 = 'positive';
                                    case 'arean'
                                          op2 = 'negative';
                              end
                              
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;                             
                              A      =  areaerp(dataux, ERP.srate, latsamp, op2) ;
                              VALUES(b,ch)  = A;
                              
                        elseif strcmpi(op1,'50arealat')  % obsolete
                              %
                              % get 50% area latency (old)
                              %
                              
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;                              
                              [aaaxxx L]   =  areaerp(dataux, ERP.srate,latsamp) ;
                              VALUES(b,ch) = ((L-1)/fs + ERP.xmin)*1000; % 50 % area latency (temporary)
                              
                        elseif strcmpi(op1,'fareatlat') || strcmpi(op1,'fninteglat') ||  strcmpi(op1,'fareaplat') || strcmpi(op1,'fareanlat')
                              
                              %
                              % get fractional area latency
                              %
                              if frac<0 || frac>1
                                    error('ERPLAB says: error at geterpvalues(). Fractional area value must be between 0 and 1')
                              end
                              switch op1
                                    case 'fareatlat'
                                          op2 = 'total';
                                    case 'fninteglat'
                                          op2 = 'integral'; % default
                                    case 'fareaplat'
                                          op2 = 'positive';
                                    case 'fareanlat'
                                          op2 = 'negative';
                              end
                              
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;
                              [aaaxxx L]   =  areaerp(dataux, ERP.srate,latsamp, op2,1, frac, fracmearep) ;
                              VALUES(b,ch) = ((L-1)/fs + ERP.xmin)*1000; % frac area latency
                              
                        elseif strcmpi(op1,'fpeaklat')
                              
                              %
                              % get fractional "peak" latency
                              %
                              if frac<0 || frac>1
                                    error('ERPLAB says: error at geterpvalues(). Fractional peak value must be between 0 and 1')
                              end
                              
                              %                               try
                              blv = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              
                              dataux  = ERP.bindata(chanArray(ch),latsamp(1)-sampeak:latsamp(2)+sampeak, binArray(b)) - blv;  % base line was substracted!
                              timex   = ERP.times(latsamp(1)-sampeak:latsamp(2)+sampeak);
                              
                              if localopt==1 %0=writes a NaN when local peak is not found.  1=export absolute peak when local peak is not found.
                                    localoptstr = 'abs';
                              else
                                    localoptstr = 'NaN';
                              end
                              if fracmearep==1 %0=writes a NaN when local peak is not found.  1=export absolute peak when local peak is not found.
                                    fracmearepstr = 'abs';
                              else
                                    fracmearepstr = 'NaN';
                              end
                              
                              [aaaxxx latpeak latfracpeak] = localpeak(dataux, timex, 'Neighborhood',sampeak, 'Peakpolarity', polpeak, 'Sfactor', 4, 'Measure','fraclat',...
                                    'Peakreplace', localoptstr, 'Fraction', frac, 'Fracpeakreplace', fracmearepstr);
                              
                              if isempty(aaaxxx)
                                    error('Peak-related measurement failed...')
                              end
                              
                              worklate{b,ch} = latpeak; % peak
                              VALUES(b,ch)   = latfracpeak; % fractional peak
                              
                        elseif strcmpi(op1,'areazt') || strcmpi(op1,'areazp') || strcmpi(op1,'areazn')
                              
                              %
                              % get area (automatic limits)
                              %
                              switch op1
                                    case 'areazt'
                                          op2 = 'autot';
                                    case 'areazp'
                                          op2 = 'autop';
                                    case 'areazn'
                                          op2 = 'auton';
                              end
                              
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;                            
                              [A L il]  =  areaerp(dataux, ERP.srate,latsamp, op2, coi) ;
                              worklate{b,ch} = ((il-1)/fs + ERP.xmin)*1000; % integration limits
                              VALUES(b,ch)   = A;                              
                        elseif strcmpi(op1,'errorbl') % for Rick Addante
                              %
                              % get standard deviation
                              %
                              if isempty(ERP.binerror)
                                    error('ERPLAB says: Rick, the data field for standard deviation is empty!')
                              end
                              
                              dataux = ERP.bindata; % temporary store for data field
                              ERP.bindata = ERP.binerror;  % error data to data
                              blv  = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value                              
                              VALUES(b,ch) = mean(ERP.bindata(chanArray(ch),latsamp(1):latsamp(2), binArray(b))) - blv;
                              ERP.bindata  = dataux; % recover original data                              
                        elseif strcmpi(op1,'ninteg')
                              %
                              % get numerical integration
                              %
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;
                              A  =  areaerp(dataux, ERP.srate,latsamp, 'integral', coi) ;
                              VALUES(b,ch)  = A;                              
                        elseif strcmpi(op1,'nintegz')
                              %
                              % get numerical integration (automatic limits)
                              %
                              blv    = blvalue(ERP, chanArray(ch), binArray(b), blc); % baseline value
                              dataux = ERP.bindata(chanArray(ch), :, binArray(b)) - blv;
                              [A L il]  =  areaerp(dataux, ERP.srate,latsamp, 'auto', coi) ;                              
                              worklate{b,ch} = ((il-1)/fs + ERP.xmin)*1000; % integratin limits
                              VALUES(b,ch)   = A;
                        end
                  end
            end
      end
catch
      serr = lasterror;
      varargout{1} = serr.message;
      varargout{2} = [];
      return
end

%
% Creates Output
%
if ~condf
      varargout{1} = VALUES;
elseif condf
      varargout{1} = VALUES;
      varargout{2} = worklate;
else
      error('ERPLAB says: error at geterpvalues.  Too many output arguments!')
end

% % %---------------------------------------------------------------------------------------------------
% % %-----------------base line mean value--------------------------------------------------------------
% % function blv = blvalue(ERP, chan, bin, blcorr)
% % 
% % %
% % % Baseline assessment
% % %
% % if ischar(blcorr)
% %       
% %       if ~strcmpi(blcorr,'no') && ~strcmpi(blcorr,'none')
% %             
% %             if strcmpi(blcorr,'pre')
% %                   bb = find(ERP.times==0);    % zero-time locked
% %                   aa = 1;
% %             elseif strcmpi(blcorr,'post')
% %                   bb = length(ERP.times);
% %                   aa = find(ERP.times==0);
% %             elseif strcmpi(blcorr,'all') || strcmpi(blcorr,'whole')
% %                   bb = length(ERP.times);
% %                   aa = 1;
% %             else
% %                   blcnum = str2num(blcorr); % in ms
% %                   
% %                   %
% %                   % Check & fix baseline range
% %                   %
% %                   if blcnum(1)<ERP.xmin*1000
% %                         blcnum(1) = ERP.xmin*1000; %ms
% %                   end
% %                   if blcnum(2)>ERP.xmax*1000
% %                         blcnum(2) = ERP.xmax*1000; %ms
% %                   end
% %                   
% %                   [xxx, cindex] = closest(ERP.times, blcnum); % 04/21/2011
% %                   aa = cindex(1); % ms to sample pos
% %                   bb = cindex(2); % ms to sample pos
% %             end
% %             blv = mean(ERP.bindata(chan,aa:bb, bin));
% %       else
% %             blv = 0;
% %       end
% % else      
% %       %
% %       % Check & fix baseline range
% %       %
% %       if blcorr(1)<ERP.xmin*1000
% %             blcorr(1) = ERP.xmin*1000; %ms
% %       end
% %       if blcorr(2)>ERP.xmax*1000
% %             blcorr(2) = ERP.xmax*1000; %ms
% %       end
% %       [xxx, cindex] = closest(ERP.times, blcorr); % 04/21/2011
% %       aa = cindex(1); % ms to sample pos
% %       bb = cindex(2); % ms to sample pos
% %       blv = mean(ERP.bindata(chan,aa:bb, bin));
% % end


