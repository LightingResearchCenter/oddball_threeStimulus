% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2011

function [msg2show msgori mcolor] = readmsg2end(textin)
try
      if nargin<1
            p = which('eegplugin_erplab');
            p = p(1:strfind(p,'eegplugin_erplab.m')-1);
            filename = fullfile(p,'functions','msg2end.txt');
            fid = fopen(filename);
            k=1; msgArray{1,1} = [];
            while ~feof(fid) && k<100
                  msg = textscan(fid, '%[^\n]',1);
                  msgArray{k,1} = char(msg{:});
                  k=k+1;
            end
            fclose(fid);
      else
            if iscell(textin)
                  msgArray = textin;
            else
                  msgArray = cellstr(textin);
            end
      end
      i=1; j=1; findcol=1; msgcell2show{1,1} = []; msgcellori{1,1} = []; mcolorx = [];
      for g=1:length(msgArray)
            msg = msgArray{g};
            auxmsg = msg;
            if isempty(msg)
                  msg = NaN;
            else
                  if findcol==1
                        rexp = '<\s*\[\s*(\d\.*\d*)*\s*(\d\.*\d*)*\s*(\d\.*\d*)*\s*]\s*>';
                        colval = regexp(msg, rexp, 'tokens');
                        if length([colval{:}])==3
                              c1 = str2double(colval{1}{1});
                              c2 = str2double(colval{1}{2});
                              c3 = str2double(colval{1}{3});
                              mcolorx = [c1 c2 c3];
                              if length(mcolorx)==3
                                    msg = Inf;
                                    findcol = 0;
                              else
                                    mcolorx = [];
                              end
                        end
                  end
            end
            if ~isnumeric(msg)
                  evalcom = regexp(msg, '\s*eval\(\''(.)*\''\)', 'tokens');
                  evalcom = strtrim(char(evalcom{:}));
                  if ~isempty(evalcom);
                        eval(evalcom);
                        msg = NaN;
                  end
            end
            if ~isnumeric(msg)
                  msgcellori{j,1} = msg;
                  try
                        msgx = eval(msg) ; % eval the message
                        msg = msgx;
                  catch
                  end
                  msgcell2show{i,1} = sprintf('%s\n', msg); %[msg '\n'];
                  i=i+1;
                  j=j+1;
            elseif isnan(msg)
                  msgcellori{j,1} = auxmsg;
                  j=j+1;
            end
      end
      if isempty([msgcell2show{:}])
            msg2show = 'COMPLETE';% default
      else
            msg2show = char(cellstr(msgcell2show));
      end
      if isempty([msgcellori{:}])
            msgori = 'COMPLETE';% default
      else
            msgori = char(cellstr(msgcellori));
      end
      if isempty(mcolorx)
            mcolor = [0 0 1]; % default blue
      else
            mcolor = mcolorx;
      end
catch
      msg2show = 'COMPLETE'; % default
      msgori   = 'COMPLETE'; % default
      mcolor   = [0 0 1];    % default blue
end

