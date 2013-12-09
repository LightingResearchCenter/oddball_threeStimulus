% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2011
function erplabamnesia(warningop)

if nargin<1
    warningop = 0;
end
if warningop
%Warning Message
    question   = ['Resetting ERPLAB''s working memory will\n'...
        'Clear all memory and cannot be recovered\n'...
        'Do you want to continue anyway?'];
    title      = 'ERPLAB: Reset ERPLAB''s working memory Confirmation';
    button     = askquest(sprintf(question), title);
    
    if ~strcmpi(button,'yes')
        disp('User selected Cancel')
        return
    end
end

p = which('eegplugin_erplab');
p = p(1:findstr(p,'eegplugin_erplab.m')-1);
recycle on;
mfile = fullfile(p,'memoryerp.erpm');
if exist(mfile, 'file')
        v=load(mfile, '-mat');
        erplabver  = v.erplabver;
        erplabrel  = v.erplabrel;  
        mshock     = v.mshock + 1;
        delete(mfile)
        fprintf('\n*** ERPLAB WARNING: ERPLAB''s working memory was wiped out. Default values will be used.\n\n')
        ColorB = [170 180 195]/255;
        ColorF = [0 0 0];
        save(fullfile(p,'memoryerp.erpm'), 'erplabrel','erplabver', 'ColorB', 'ColorF', 'mshock'); % saves erplab version
        if mshock>=10 && rand>0.8
                fprintf('\n\nIs it not enough???\n\n')
        end
else
        fprintf('\n*** ERPLAB WARNING: ERPLAB''s working memory file does not existe.\n')
        fprintf('\n*** \t\tIt is strongly recommended you restart EEGLAB. \n')
end
recycle off