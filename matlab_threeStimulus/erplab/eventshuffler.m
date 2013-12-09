% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009
%
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

function EEG = eventshuffler(EEG, valueatfield, specfield)
if nargin<1
        help eventshuffler
        return
end
if nargin<2
        specfield = 0; % codes {default}
end
if specfield==0
        %
else
        if isfield(EEG, 'EVENTLIST')
                if isempty(EEG.EVENTLIST)
                        %fprintf('\nWARNING: eventshuffler() did not shuffle your bins.\n');
                        error('eventshuffler() cannot work on an empty EEG.EVENTLIST structure.');
                end
        else
                %fprintf('\nWARNING: eventshuffler() did not shuffle your bins.\n');
                error('eventshuffler() cannot shuffle bins without EEG.EVENTLIST structure.');
        end
end
if isnumeric(valueatfield) && length(valueatfield)>1
        if specfield==0
                if isfield(EEG, 'EVENTLIST') && ~isempty(EEG.EVENTLIST)
                        indx = find(ismember([EEG.EVENTLIST.eventinfo.code], valueatfield));
                else
                        indx = find(ismember([EEG.event.type], valueatfield));
                end
        else
                k=1;
                for a=1:length(EEG.EVENTLIST.eventinfo)
                        if nnz(ismember( valueatfield, [EEG.EVENTLIST.eventinfo(a).bini]))>0;
                                indx(k) = a;
                                k = k+1;
                        end
                end
        end
elseif ischar(valueatfield) && strcmpi(valueatfield, 'all')
        if specfield==0
                indx = 1:length(EEG.EVENTLIST.eventinfo);
        else
                k=1;
                for a=1:length(EEG.EVENTLIST.eventinfo)
                        if nnz([EEG.EVENTLIST.eventinfo(a).bini]>0)>0
                                indx(k) = a;
                                k = k+1;
                        end
                end
        end
else
        %fprintf('\nWARNING: eventshuffler() did not shuffle your %ss.\n', w);
        error('ERROR:eventshuffler', ['eventshuffler() needs two or more codes to do the job.\n'...
              'For shuffling all your codes please specify ''all'' as the second input']);
end
if specfield==0 && isempty(EEG.epoch) % shuffle codes (only continuous)
        if isfield(EEG, 'EVENTLIST') && ~isempty(EEG.EVENTLIST)
                                
                
                codes  = [EEG.EVENTLIST.eventinfo(indx).code];
                ncodes = length(codes);
                p      = randperm(ncodes);
                codes  = num2cell(codes(p));
                [EEG.EVENTLIST.eventinfo(indx).code] = codes{:};
                EEG = pop_overwritevent( EEG, 'code');
        else
                codes  = [EEG.event(indx).type];
                ncodes = length(codes);
                p      = randperm(ncodes);
                codes  = num2cell(codes(p));
                [EEG.event(indx).type] = codes{:};
        end
elseif specfield==1 % shuffle bins (continuous & epoched)
        if isempty(EEG.epoch) % continuous
                nindx = length(indx);
                p     = randperm(nindx);
                rindx = indx(p);
                
                for a=1:length(indx)
                        bini   = EEG.EVENTLIST.eventinfo(rindx(a)).bini;
                        flag   = EEG.EVENTLIST.eventinfo(rindx(a)).flag;
                        enable = EEG.EVENTLIST.eventinfo(rindx(a)).enable;
                        
                        EEG.EVENTLIST.eventinfo(indx(a)).bini   = bini;
                        EEG.EVENTLIST.eventinfo(indx(a)).flag   = flag;
                        EEG.EVENTLIST.eventinfo(indx(a)).enable = enable;
                        
                        % Creates BIN LABELS
                        auxname = num2str(bini);
                        bname   = regexprep(auxname, '\s+', ',', 'ignorecase'); % inserts a comma instead blank space
                        
                        if strcmp(EEG.EVENTLIST.eventinfo(rindx(a)).codelabel,'""')
                                binName = ['B' bname '(' num2str(EEG.EVENTLIST.eventinfo(rindx(a)).code) ')']; %B#(code)
                        else
                                binName = ['B' bname '(' EEG.EVENTLIST.eventinfo(rindx(a)).codelabel ')']; %B#(codelabel)
                        end
                        
                        EEG.EVENTLIST.eventinfo(indx(a)).binlabel = binName;
                end
                
                EEG = pop_overwritevent( EEG, 'binlabel');
        else
                %fprintf('\nWARNING: eventshuffler() did not shuffle your %ss.\n', w);
                error('eventshuffler() can shuffle events in continuous data only.');                
                %nepoch = EEG.trials;
                %for e=1:nepoch
                %
                %end
        end
else
        %fprintf('\nWARNING: eventshuffler() did not shuffle your %ss.\n', w);
        error('eventshuffler() can shuffle event types in continuous data only.');
end

% EEG = eeg_checkset( EEG );



