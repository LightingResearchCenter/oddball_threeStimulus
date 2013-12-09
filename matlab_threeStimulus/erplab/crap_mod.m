% crap.m (alpha version). Reject commonly recorded artifactual potentials (c.r.a.p.)
%
% There are a number of common artifacts that you will see in nearly every EEG data file. These
% include eyeblinks, slow voltage changes (caused mostly by skin potentials), muscle activity
% (from moving the head or tensing up the muscles in the face or neck), horizontal eye movements,
% and various types of C.R.A.P. (Commonly Recorded Artifactual Potentials).
%
% Although we usually perform artifact rejection on the segmented data, it's a good idea to
% examine the raw unsegmented EEG data first. You can usually identify patterns of artifacts,
% make sure there were no errors in the file, etc., more easily with the raw data [1].
%
% crap.m allows you to automatically identify large peak-to-peak differences or extreme amplitude
% values, within a moving window, across your continuous EEG dataset. After performing crap.m,
% artifactual segments will be rejected and replaced by a 'boundary' event code.
%
% USAGE
%
% EEG = crap(EEG, ampth, winms, stepms,chanArray, option);
%
% Input:
%
% EEG         - continuous EEG dataset (EEGLAB's EEG structure)
% ampth       - 1 single value for peak-to-peak threshold within the moving window
%             - 2 values for extreme thresholds within the moving window, e.g [-200 200] or [-150 220]
% winms       - moving window width in msec (default 2000 ms)
% stepms      - moving window step (default 1000 ms)
% chanArray   - channel index(ices) to look for c.r.a.p.  (default all channels)
%
% Option is related to filtering:
% 'filter'    - add 2 values indicting the frequency cutoff. For instance, [8 12]
%               Also, you may add a string that will call a pre-defined
%               cutoff. For instamce: 'rdc','rdelta','rtheta','ralpha','rbeta','rgamma'
%               or 'dc','delta','theta','alpha','beta','gamma'
%
%              'rdc'     - will remove the DC offset. Equivalent to [0 0] (remove the mean value). Default.
%              'rdelta'  - will remove delta band. Equivalent to [4 0.1]
%              'rtheta'  - will remove theta band. Equivalent to [8 4]
%              'ralpha'  - will remove alpha band. Equivalent to [13 8]
%              'rbeta'   - will remove beta band. Equivalent to [30 13]
%              'rgamma'  - will remove gamma band. Equivalent to [80 30]
%
%              'dc'     - will remove everything but the DC offset. Equivalent to [Inf Inf] (only keep the mean value).
%              'delta'  - will remove everything but  delta band. Equivalent to [0.1 4]
%              'theta'  - will remove everything but  theta band. Equivalent to [4 8]
%              'alpha'  - will remove everything but  alpha band. Equivalent to [8 13]
%              'beta'   - will remove everything but  beta band.  Equivalent to [13 30]
%              'gamma'  - will remove everything but  gamma band. Equivalent to [30 80]
%
%              'off' or 'no' -  will disable any filter process.
%
% 'forder'    - 1 value indicting the order (number of points) of the FIR filter to be used.
%               When 'filter' is specified (other than 'rdc' or 'dc') the default value is 26.
%
%
% Output:
%
% EEG         - continuous EEG dataset without c.r.a.p. (EEGLAB's EEG structure)
%
%
% Example 1:
% Reject segment of data where the amplitude (mean value removed) is >= than 200 uV or <= than -300 uV.
% Use a moving window of 2000 ms, 1000 ms step, exploring channels 1 to 12.
%
% >> EEG = crap(EEG, [-300 200], 2000, 1000, 1:12);
%
% Example 2:
% Reject segment of data where the instantaneous amplitude is >= than 120 uV or <= than -120 uV.
% Use a moving window of 2000 ms, 1000 ms step, exploring channels 1 to 12.
%
% >> EEG = crap(EEG, [-120 120], 2000, 1000, 1:12, 'filter','off');
%
% Example 3:
% Reject segment of data where the peak-to-peak amplitude is >= 240 uV. Note that mean values does not matter in this case.
% Use a moving window of 2000 ms, 1000 ms step, exploring channels 1 to 12.
%
% >> EEG = crap(EEG, 240, 2000, 1000, 1:12);
%
% Example 4:
% Reject segment of data where beta peak-to-peak amplitude is higher than 10 uV.
% Use a moving window of 1000 ms, 500 ms step, exploring channels 40 and 42.
%
% >> EEG = crap(EEG, 10, 1000, 500, [40 42], 'filter',[13 30]);
%    or
% >> EEG = crap(EEG, 10, 1000, 500, [40 42], 'filter','beta');
%
% Example 5:
% Reject segment of data where beta peak-to-peak amplitude is higher than 10 uV.
% Use a moving window of 1000 ms, 500 ms step, exploring channels 40 and 42.
%
% >> EEG = crap(EEG, 10, 1000, 500, [40 42], 'filter',[13 30]);
%    or
% >> EEG = crap(EEG, 10, 1000, 500, [40 42], 'filter','beta');
%
% Example 6:
% Reject segment of data where any peak-to-peak activity, but alpha, is higher than 80 uV.
% Use a moving window of 3000 ms, 1000 ms step, exploring channels all 60 channels.
% Use a 60th order filter.
%
% >> EEG = crap(EEG, 80, 3000, 100, 1:60, 'filter',[13 8], 'forder', 60);
%    or
% >> EEG = crap(EEG, 80, 3000, 100, 1:60, 'filter','ralpha', 'forder', 60);
%
% Example 7 (rare):
% Reject segment of data where the DC offset is >= than 15000 or <= than -15000 uV.
% Use a moving window of 2000 ms, 1000 ms step, exploring channels 1 to 10 and 16 to 32.
%
% >> EEG = crap(EEG, [-15000 15000], 2000, 1000, [1:10 16:32], 'filter','dc');
%
%
%
% Reference:
% [1] ERP Boot Camp: Data Analysis Tutorials. Emily S. Kappenman, Marissa L. Gamble, and Steven J. Luck. UC Davis
%
%
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009

%b8d3721ed219e65100184c6b95db209bb8d3721ed219e65100184c6b95db209b
%
% ERPLAB Toolbox
% Copyright � 2007 The Regents of the University of California
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
%
%
% Debugging :
% bug at ampth with 2 values was fixed. JLC & J.Kreither
% remove dc option was added
%
function [NaN_indices, vDiffOut] = crap_mod(EEG, ampth, winms, stepms, chanArray, varargin)

if nargin<1
      help crap
      return
end
if nargin<2
      error('ERPLAB says: crap.m needs at least 2 inputs.')
end
if nargin<5
      chanArray = 1:EEG.nbchan; % channel to look for c.r.a.p.
end
if nargin<4
      stepms = 1000; % ms ; moving window step
end
if nargin<3
      winms = 2000; % ms ; moving window width
end
if length(EEG)>1
      msgboxText =  'Unfortunately, this function does not work with multiple datasets';
      error(msgboxText)
end
if ~isempty(EEG.epoch)
      msgboxText =  'crap() only works for continuous datasets.';
      error(msgboxText)
end
if numel(ampth)~=1 && numel(ampth)~=2
      error('ERPLAB says: for threshold amplitude you must specify 1 or 2 values only.')
end

p = inputParser;
p.FunctionName  = mfilename;
p.CaseSensitive = false;
p.addRequired('EEG', @isstruct);
p.addRequired('ampth', @isnumeric);
p.addRequired('winms', @isnumeric);
p.addRequired('stepms', @isnumeric);
p.addRequired('chanArray', @isnumeric);
p.addParamValue('filter', 'dc');
p.addParamValue('forder', 26, @isnumeric);

p.parse(EEG, ampth, winms, stepms, chanArray, varargin{:});

fcutoff = p.Results.filter;
forder  = p.Results.forder;
if mod(forder,2)~=0
      forder = forder+1;
      fprintf('filter order was changed to an even number (%g) because of the forward-reverse filtering.\n', forder)
end
if ischar(fcutoff)
      %
      % Pre-defined bands
      %
      % Removing
      if strcmpi(fcutoff,'rdc') || strcmpi(fcutoff,'rmean')
            fcutoff = [0 0];           % remove mean
      elseif strcmpi(fcutoff,'rdelta') % remove delta
            fcutoff = [4 0.1];
      elseif strcmpi(fcutoff,'rtheta') % remove theta
            fcutoff = [8 4];
      elseif strcmpi(fcutoff,'ralpha') % remove alpha
            fcutoff = [13 8];
      elseif strcmpi(fcutoff,'rbeta')  % remove beta
            fcutoff = [30 13];
      elseif strcmpi(fcutoff,'rgamma') % remove gamma
            fcutoff = [80 30];
            % keeping
      elseif strcmpi(fcutoff,'dc') || strcmpi(fcutoff,'mean') % keep only the mean
            fcutoff = [inf inf];
      elseif strcmpi(fcutoff,'delta') % keep only delta
            fcutoff = [0.1 4];
      elseif strcmpi(fcutoff,'theta') % keep only theta
            fcutoff = [4 8];
      elseif strcmpi(fcutoff,'alpha') % keep only alpha
            fcutoff = [8 13];
      elseif strcmpi(fcutoff,'beta')  % keep only beta
            fcutoff = [13 30];
      elseif strcmpi(fcutoff,'gamma') % keep only gamma
            fcutoff = [30 80];
      elseif strcmpi(fcutoff,'off') || strcmpi(fcutoff,'no') % turn off filtering
            fcutoff = [];
      else
            error('error: unknown pre-defined band. Try with a numeric range of frequencies.')
      end
else
      if ~isempty(fcutoff)
            if numel(fcutoff)~=2
                  error('error: 2 values are needed for cutoff.')
            end
            if (~isinf(fcutoff(1)) && ~isinf(fcutoff(2))) ||...
                        (~isnan(fcutoff(1)) && ~isnan(fcutoff(2)))
                  if fcutoff(1)~=0 && ~isinf(fcutoff(1)) && fcutoff(1)==fcutoff(2)
                        error('error: cutoff must have 2 different values. Except for [0 0] or [Inf Inf]. See help crap')
                  end
                  if fcutoff(1)<0 || fcutoff(2)<0
                        error('error: Cutoff values must be >= 0, or Inf')
                  end
                  if fcutoff(1)>EEG.srate/2 || fcutoff(2)>EEG.srate/2
                        error('error: frequency cutoff is higher than nyquist frequency (EEG.srate/2)')
                  end
            else
                  if (isinf(fcutoff(1)) && ~isinf(fcutoff(2))) ||...
                              (~isinf(fcutoff(1)) && isinf(fcutoff(2))) ||...
                              (isnan(fcutoff(1)) || isnan(fcutoff(2)))
                        error('error: cutoff can not have NaN or Inf number. Only [Inf Inf] is allowed for evaluating the mean value.')
                  end                  
            end
      else
            error('error: cutoff can not be empty')
      end
end
% disp('Working...')
nchan = length(chanArray);
fs    = EEG.srate;
dursam1 = EEG.pnts;

%
% for searching boundaries inside EEG.event.type
%
if length(EEG.event)<1
      fprintf('\ncrap.m.m did not found remaining event codes.\n')
      return
end
if ischar(EEG.event(1).type)
      codebound = {EEG.event.type}; %strings
      indxbound = strmatch('boundary', codebound, 'exact');
else
      indxbound = [];
end
if ~isempty(indxbound)
      timerange = [ EEG.xmin*1000 EEG.xmax*1000 ];
      if timerange(1)/1000~=EEG.xmin || timerange(2)/1000~=EEG.xmax
            posi = round( (timerange(1)/1000-EEG.xmin)*EEG.srate )+1;
            posf = min(round( (timerange(2)/1000-EEG.xmin)*EEG.srate )+1, EEG.pnts );
            pntrange = posi:posf;
      end
      if exist('pntrange')
            latebound = [ EEG.event(indxbound).latency ] - 0.5 - pntrange(1) + 1;
            latebound(find(latebound>=pntrange(end) - pntrange(1))) = [];
            latebound(find(latebound<1)) = [];
            latebound = [0 latebound pntrange(end) - pntrange(1)];
      else
            latebound = [0 [ EEG.event(indxbound).latency ] - 0.5 EEG.pnts ];
      end
      latebound = round(latebound);
else
      latebound = [0 EEG.pnts];
      % fprintf('\nWARNING: boundary events were not found.\n');
      % fprintf('WARNING: crap.m will be applied over the full range of data.\n\n');
end
nibound   = length(latebound);
winpnts  = floor(winms*fs/1000); % to samples
stepnts  = floor(stepms*fs/1000);% to samples
winrej = []; k=1;

for ch=1:nchan
      q=1;
      while q<=nibound-1  % segments among boundaries
            bp1   = latebound(q)+1;
            bp2   = latebound(q+1);
            % fprintf('Exploring channel %g, segment %g to %g (in samples)...\n', chanArray(ch), bp1, bp2);
            datax = EEG.data(chanArray(ch), bp1:bp2);
            if ~isempty(fcutoff)
                  if fcutoff(1)~=fcutoff(2)
                        if length(bp1:bp2)>3*forder
                              % FIR coefficients
                              [b, a, labelf] = filter_tf(1, forder, fcutoff(2), fcutoff(1), EEG.srate);
                              if fcutoff(1)>0 % option to remove dc value is only for high-pass filtering
                                    %
                                    % Removes DC
                                    %
                                    % fprintf('Removing DC bias before high/band pass filtering...\n');
                                    datax = detrend(datax, 'constant');
                              end
                              fprintf('%s filtering data; cutoff = [%.1f %.1f]...\n', labelf, fcutoff(1), fcutoff(2));                              
                              % FIR lowpass
                              % FIR highpass
                              % FIR bandpass
                              % FIR notch
                              if isdoublep(datax)
                                    datax = filtfilt(b,a, datax);
                              else
                                    datax = single(filtfilt(b,a, double(datax)));
                              end
                        else
                              fprintf('\nWARNING: EEG segment from sample %d to %d was not filtered\n', bp1,bp2);
                              fprintf('because number of samples must be >= 3 x filter''s order.\n\n');
                        end
                        fproblems = nnz(isnan(datax));
                        if fproblems>0
                              msgboxText = ['Oops! filter is not working properly. Data have undefined numerical results.\n'...
                                    'We strongly recommend that you change some filter parameters,\n'...
                                    'for instance, decrease filter order.'];
                              msgboxText = sprintf(msgboxText);
                              error(msgboxText);
                        end                        
                  end
            end
            
            %
            % Moving window
            %
            endpoint = bp2;
            j=bp1;
            while j<=endpoint-(winpnts-1)
                  t1   = j+1;
                  t2   = j+winpnts-1;
                  datax2 = datax(t1-bp1:t2-bp1);
                  if ~isinf(fcutoff(1))
                        if (fcutoff(1)+fcutoff(2))==0
                              if j==1
                                    % fprintf('Removing DC offset from windows within segment %g, channel %g \n', q, chanArray(ch))
                              end
                              datax2 = detrend(datax2, 'constant');
                        end
                  else
                        if j==1
                              % fprintf('Isolating DC offset (mean) from windows within segment %g, channel %g \n', q, chanArray(ch))
                        end
                        datax2 = datax2-detrend(datax2, 'constant');
                  end
                  vmin = min(datax2); vmax = max(datax2);
                  if length(ampth)==1
                        vdiff   = abs(vmax - vmin);
                        if vdiff>ampth
                              winrej(k,:) = [t1 t2]; % start and end samples for rejection
                              k=k+1;
                        end
                        vDiffOut(j) = vdiff;
                  else
                        if vmin<=ampth(1) || vmax>=ampth(2)
                              winrej(k,:) = [t1 t2]; % start and end samples for rejection
                              k=k+1;
                        end
                        vDiffOut(j) = max([abs(vmin) abs(vmax)]);
                  end
                  j=j+stepnts;
            end
            q = q + 1; % next segment
      end
end
if isempty(winrej)
      % fprintf('\nCriterion was not found. No rejection was performed.\n');
      NaN_indices = logical(zeros(length(EEG.data)));
else
      % Selects not overlapping start and end samples
      winrej = sort(winrej,1);
      [aa1 bb1] = unique( winrej(:,1),'first');
      [aa2 bb2] = unique( winrej(:,2),'last');
      winrej = [winrej(bb1,1) winrej(bb2,2)];
      a = winrej(1,1); winrej2(1,:) = [a winrej(end,2)]; m=1;
      for j=2:size(winrej,1)
            if abs(winrej(j,2)-winrej(j-1,2))>winpnts
                  b = winrej(j-1,2);
                  winrej2(m,:) = [a b];
                  a = winrej(j,1);
                  m=m+1;
            end
      end
      if winrej2(end,2)~=winrej(end,2)
            winrej2(m,:) = [a winrej(end,2)];
      end
      
      % rejects      
      NaN_indices = zeros(length(EEG.data),1);
      for rejEp = 1 : size(winrej2,1)
          NaN_indices(winrej(rejEp,1):winrej(rejEp,2)) = 1;          
      end
      NaN_indices = logical(NaN_indices);
      
      
      
      
      % PETTERI: The following checks are removed, and we only need the
      % indices found by the routine
      
      %{
      
      EEG = eeg_eegrej(EEG, winrej2);
      
      EEG = eeg_checkset( EEG );
      if length(EEG.event)>=1
            if EEG.event(end).latency>EEG.pnts
                  EEG = pop_editeventvals(EEG,'delete',length(EEG.event));
                  EEG = eeg_checkset( EEG );
            end
      end
      if length(EEG.event)>=1
            if EEG.event(1).latency<1
                  EEG = pop_editeventvals(EEG,'delete',1);
                  EEG = eeg_checkset( EEG );
            end
      end
      
      EEG = delshortseg(EEG,'boundary',100); % segments shorter than 100 msec will be deleted
      dursam2 = EEG.pnts;
      fprintf([repmat('-',1,60) '\n']);
      fprintf('Cost:\nYour dataset was shortened %.1f percent.\n', 100-100*(dursam2/dursam1))
      fprintf([repmat('-',1,60) '\n\n']);
      %}
end