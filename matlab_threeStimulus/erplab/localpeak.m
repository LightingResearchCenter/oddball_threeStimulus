% Author: Javier Lopez-Calderon & Steve Luck
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2011

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

function [measure, latlocalpeak, latfracpeak, errorcode] = localpeak(datax, timex, varargin)

measure      = [];
vlocalpf     = [];
latlocalpeak = [];
latfracpeak  = [];
errorcode    = [];

%
% Parsing inputs
%
p = inputParser;
p.FunctionName  = mfilename;
p.CaseSensitive = false;
p.addRequired('datax', @isnumeric);
p.addRequired('timex', @isnumeric);

p.addParamValue('Neighborhood', 1); % erpset index or input file
p.addParamValue('Peakpolarity', 1, @isnumeric);
p.addParamValue('Multipeak', 'off', @ischar);
p.addParamValue('Fraction', [], @isnumeric);
p.addParamValue('Sfactor', 1, @isnumeric);
p.addParamValue('Measure', 'amplitude', @ischar); % 'amplitude', 'peaklat', 'fraclat'
p.addParamValue('Peakreplace', 'NaN', @ischar); % alterantive amplitude (when no local peak was found)
p.addParamValue('Fracpeakreplace', 'NaN', @ischar); % alterantive amplitude (when no local peak was found)
% p.addParamValue('Altlat', 'abs', @ischar); % alterantive latency (when no local peak was found)

try
      p.parse(datax, timex, varargin{:});
catch
      errorcode = 1;
      return
end

%
% Tests flat data
%
du = unique(datax);
if length(du)<=round(length(datax)*0.1) % unique values are less than 10% of total sample values
      errorcode = 6;
      return
end
if strcmpi(p.Results.Multipeak,'on') || strcmpi(p.Results.Multipeak,'yes')
      multi = 1;
else
      multi = 0;
end

frac     = p.Results.Fraction;             % Fractional peak
npoints  = round(p.Results.Neighborhood);  % sample(s) at one side of the peak
nsamples = length(datax);
peakpol  = p.Results.Peakpolarity;
Fk       = round(p.Results.Sfactor);  % oversampling factor

%
% Oversampling
%
if Fk>1
      p1  = timex(1);
      p2  = timex(end);
      timex2   = linspace(p1,p2,Fk*(p2-p1+1));
      datax    = spline(timex, datax, timex2); % over sampled data
      timex    = timex2;
      npoints  = npoints*Fk;    % scale neighborhood
      nsamples = length(datax); % new size
end
if nsamples<=(2*npoints)
      errorcode = 1; % error. few samples.
      return
end

%
% % Starting points, npoints=neighbors, nsamples=#of samples of the window
% of interest + 2*npoints
%
a = npoints + 1;
b = nsamples - npoints;

%
% Absolute peaks
%
try
      if peakpol==1 % positive peak -> finds maximum
            [vabspf posabspf] = max(datax(a:b));
      else  % negative peak -> finds minimum
            [vabspf posabspf] = min(datax(a:b));
      end
      posabspf = posabspf + npoints;
      %latabspeak = timex(posabspf); % latency for absolute peak
catch
      errorcode = 1;
      return
end

%
% local peaks
%
k = 1;
valmax = [];

if npoints>0 % LOCAL
      while a<=b
            
            avgneighborLeft  = mean(datax(a-npoints:a-1));
            avgneighborRight = mean(datax(a+1:(a+npoints)));
            prea  = datax(a-1);
            posta = datax(a+1);
            
            if peakpol==1 % maximum for positives
                  if datax(a)>avgneighborLeft && datax(a)>avgneighborRight && datax(a)>prea && datax(a)>posta
                        valmax(k) = datax(a);
                        posmax(k) = a;
                        k=k+1;
                  end
            else  % minimum for negatives
                  if datax(a)<avgneighborLeft && datax(a)<avgneighborRight && datax(a)<prea && datax(a)<posta
                        valmax(k) = datax(a);
                        posmax(k) = a;
                        k=k+1;
                  end
            end
            a = a+1;
      end
      
      if ~isempty(valmax)
            
            if length(unique(valmax))==1 && length(valmax)>1 % this is when more than one sample meets the criterias for a local peak (e.g. saturated segments)
                  poslocalpf   = round(median(posmax));    % position of local peak
                  %latlocalpeak = timex(poslocalpf);        % latency for absolute peak
                  vlocalpf     = unique(valmax);           % value of local peak
                  
            elseif length(unique(valmax))>1
                  
                  %if multi
                  %        vlocalpf = valmax; % values for multiple local peaks
                  %        [tfxxclc, poslocalpf] = ismember(valmax, datax);
                  %              latlocalpeak    = timex(poslocalpf);   % latencies for multiple local peaks
                  %        else
                  %        if peakpol==1 % positive peak -> finds maximum
                  %                [vlocalpf indx] = max(valmax);
                  %        else  % negative peak -> finds minimum
                  %                [vlocalpf indx] = min(valmax);
                  %        end
                  %        poslocalpf      = posmax(indx);
                  %        latlocalpeak    = timex(poslocalpf);   % latency for local peak
                  %end
                  
                  if multi
                        vlocalpf = valmax;
                        [tf poslocalpf] = ismember(valmax, array);
                        
                        
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                  else
                        
                        if peakpol==1
                              [vlocalpf indx] = max(valmax);
                              poslocalpf      = posmax(indx);
                        else
                              [vlocalpf indx] = min(valmax);
                              poslocalpf      = posmax(indx);
                        end
                  end
                  
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%CHANGED%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  
                  
                  
            else
                  vlocalpf     = valmax;
                  poslocalpf   = posmax;
                  %latlocalpeak = timex(poslocalpf);   % latency for local peak
            end
            
            latlocalpeak = timex(poslocalpf);   % latency for local peak
            ltypeoutput  = 1; % local peak
      else
            if strcmpi(p.Results.Peakreplace,'abs') % replaces with abs peak values
                  vlocalpf     = vabspf;
                  poslocalpf   = posabspf;
                  latlocalpeak = timex(posabspf);
                  
                  ltypeoutput = 2; % abs peak
            else % replace with NaN
                  vlocalpf     = NaN;
                  poslocalpf   = NaN;
                  latlocalpeak = NaN;
                  
                  ltypeoutput = 3; % NaN instead of a peak
            end
      end
      
else % if no neighbors then ABSOLUTE peak is taken.
      vlocalpf     = vabspf;
      poslocalpf   = posabspf;
      latlocalpeak = timex(posabspf);
      ltypeoutput = 2; % abs peak
end

%
% Fractional peak latency
%
if strcmpi(p.Results.Measure,'fraclat') % fractional latency assessment
      if ~isempty(frac)
            if frac>0 && ltypeoutput~=3
                  
                  %if ~isnan(poslocalpf)
                  a = poslocalpf; % this might be local or absolute...
                  while a>0
                        currval  = datax(a);
                        if (peakpol==1 && currval<=vlocalpf*frac) || (peakpol==0 && currval>=vlocalpf*frac) % maximum for positives; miniumum for negative peak
                              posfrac = a;
                              latfracpeak = timex(posfrac);   % latency for fractional peak
                              break
                        end
                        a = a-1;
                  end
                  if isempty(latfracpeak) %&& strcmpi(p.Results.Fracpeakreplace,'abs') % replaces with frac abs peak values
                        latfracpeak = NaN;
                  end
                  
                  %else
                  %        latfracpeak = NaN;
                  %end
            elseif frac==0 && ltypeoutput~=3
                  %posfrac = 1;
                  latfracpeak = timex(1);   % latency for fractional peak
            else
                  latfracpeak = NaN;   % latency for fractional peak
                  
            end
      else
            errorcode = 1; % no fractional value was specified
            return
      end
end

%
% Output(s)
%
if nargout==1
      if strcmpi(p.Results.Measure,'amplitude')
            measure = vlocalpf;
      elseif strcmpi(p.Results.Measure,'peaklat')
            measure = latlocalpeak;
      elseif strcmpi(p.Results.Measure,'fraclat')
            measure = latfracpeak;
      end
else
      measure = vlocalpf;
end