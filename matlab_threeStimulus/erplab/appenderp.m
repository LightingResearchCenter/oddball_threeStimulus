% - This function is part of ERPLAB Toolbox -
%

function [ERP serror] = appenderp(ALLERP,indx, prefixes)

ERP = [];
serror = 0;

if nargin<3
      prefixes=[];
end
if nargin<2
      indx=1:length(ALLERP);
end
if nargin<1
      error('ERPLAB says: appenderp.m needs at least 1 input.')
end

erpname   = [];
filename  = [];
filepath  = [];
workfiles = [];
subject   = '';
bindata   = [];
binerror  = [];
nbin      = 0;
accepted  = [];
rejected  = [];
invalid   = [];
arflags   = [];
bindescr  = [];
history   = '';
nerp    = length(indx);

for i=1:nerp      
      workfiles = [workfiles ALLERP(indx(i)).workfiles];
      subject   = [subject ALLERP(indx(i)).subject];
      
      if i==1
            nchan      = ALLERP(indx(i)).nchan;
            pnts       = ALLERP(indx(i)).pnts;
            srate      = ALLERP(indx(i)).srate;
            xmin       = ALLERP(indx(i)).xmin;
            xmax       = ALLERP(indx(i)).xmax;
            times      = ALLERP(indx(i)).times;
            chanlocs   = ALLERP(indx(i)).chanlocs;
            ref        = ALLERP(indx(i)).ref;
            bindata    = ALLERP(indx(i)).bindata;
            binerror   = ALLERP(indx(i)).binerror;
      else
            sra = size(bindata,1);
            srb = size(ALLERP(indx(i)).bindata,1);
            sca = size(bindata,2);
            scb = size(ALLERP(indx(i)).bindata,2);
            
            if sra~=srb
                  serror=2; % channel size is diff
                  return
            end
            if sca~=scb
                  serror=3; % channel size is diff
                  return
            end
            
            bindata  = cat(3, bindata, ALLERP(indx(i)).bindata);
            binerror = cat(3, binerror, ALLERP(indx(i)).binerror);
      end
      
      nbin     = nbin +  ALLERP(indx(i)).nbin;
      accepted = [accepted ALLERP(indx(i)).ntrials.accepted];
      rejected = [rejected ALLERP(indx(i)).ntrials.rejected];
      invalid  = [invalid ALLERP(indx(i)).ntrials.invalid];
      arflags  = cat(1,arflags, ALLERP(indx(i)).ntrials.arflags);
      
      if isempty(prefixes)
            bindescr = [bindescr ALLERP(indx(i)).bindescr];
      else
            auxdescr  = char(ALLERP(indx(i)).bindescr');
            auxprefix = repmat([prefixes{i} ' : '], ALLERP(indx(i)).nbin,1);
            newdescr  = cellstr(cat(2,auxprefix,auxdescr))';
            bindescr  = [bindescr newdescr];
      end
end

ERP.erpname    = erpname;
ERP.filename   = filename;
ERP.filepath   = filepath;
ERP.workfiles  = workfiles;
ERP.subject    = subject;
ERP.nchan      = nchan;
ERP.nbin       = nbin;
ERP.pnts       = pnts;
ERP.srate      = srate;
ERP.xmin       = xmin;
ERP.xmax       = xmax;
ERP.times      = times;
ERP.bindata    = bindata;
ERP.binerror   = binerror;
ERP.chanlocs   = chanlocs;
ERP.ref        = ref;
ERP.bindescr   = bindescr;
ERP.ntrials.accepted  = accepted;
ERP.ntrials.rejected  = rejected;
ERP.ntrials.invalid   = invalid;
ERP.ntrials.arflags   = arflags;
ERP.history    = history;
ERP.saved      = 'no';
ERP.isfilt     = 0;   % 1= avg was filtered or smoothed
ERP.version    = geterplabversion;

[ERP serror] = sorterpstruct(ERP);