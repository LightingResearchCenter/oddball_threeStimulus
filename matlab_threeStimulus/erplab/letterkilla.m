% Remove white space from EEG alphanumeric event codes
% Deletes non-digit character from alphanumeric even codes.
% Converts remaining codes into numeric codes.
% Unconvertibles event codes (non digit info at all) will be renamed as -88
% 
% Usage
% 
% EEG = letterkilla(EEG)
%
% Input:
% EEG     - continous dataset with alphanumeric event codes
% 
% Output
% EEG     - continous dataset with numeric event codes
%
%
%
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% January 25th, 2011
%
% Thanks to Erik St Louis and Lucas Dueffert for their valuable feedbacks.
%

function EEG = letterkilla(EEG)

if nargin<1
        help letterkilla
        return
end
%if isempty(EEG.data)
%        msgboxText = 'letterkilla() cannot read an empty dataset!';
%        title = 'ERPLAB: letterkilla() error';
%        errordlg(msgboxText,title);
%        return
%end

nevent = length(EEG.event);

if nevent<1
        msgboxText = 'Event codes were not found!';
        title = 'ERPLAB: letterkilla() error';
        errordlg(msgboxText,title);
        return
end
% if ~ischar(EEG.event(1).type)
%         %msgboxText = 'Event codes are numeric. So you do not need to run this tool.';
%         %title = 'ERPLAB: letterkilla() WARNING';
%         %errordlg(msgboxText,title);
%         fprintf('\nletterkilla() did not found alphanumeric event codes. So, events were left as they were.\n\n');
%         return
% end
fprintf('\nletterkilla() is working...\n');
EEG = wspacekiller(EEG);
for i=1:nevent
      codeaux = EEG.event(i).type;
      
      if ischar(codeaux) && ~strcmpi(codeaux, 'boundary')           
            code    = regexprep(codeaux,'\D*','', 'ignorecase'); % deletes any non-digit character
            if isempty(code)
                  code = -88;
            else
                  code = str2num(code);
            end
      else
            code = codeaux;
      end
      EEG.event(i).type =  code;
end

fprintf('letterkilla() got rid of alphabetic characters from your alphanumeric event codes.\n');
fprintf('NOTE: Event codes without any digit character, except ''boundary'', are renamed as -88 (numeric).\n\n');