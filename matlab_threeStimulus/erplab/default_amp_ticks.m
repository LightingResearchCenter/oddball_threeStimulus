% Creates default X ticks for ERPLAB's plotting GUI
%
%
% Author: Javier Lopez-Calderon & Steven Luck
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2012

function [def miny maxy] = default_amp_ticks(ERP, yrange)

def   = '-1 1';
if nargin<2 || nargout>1
      nbin  = ERP.nbin;
      ymin = zeros(1,nbin);
      ymax = ymin;
      
      for k=1:nbin
            ymin(k) = min(min(ERP.bindata(:,:,k)'));
            ymax(k) = max(max(ERP.bindata(:,:,k)'));
      end
      
      miny = min(ymin);
      maxy = max(ymax);
end
if nargin<2
      yys1 = miny*1.2;
      yys2 = maxy*1.1;      
      %[yymin*1.2 yymax*1.1]; from auto_Y at ploterps.m
else
      yys1 = yrange(1);
      yys2 = yrange(2);
end

% xxs1       = ceil(1000*ERP.xmin);
% xxs2       = floor(1000*ERP.xmax);
yystick1   = (round(yys1/10)+ 0.1)*10;
yystick2   = (round(yys2/10)+ 0.1)*10;
goags = 1; stepy = 2; L1=7; L2=15;
w=1;
while goags && w<=100
      %stepy
      ytickarray = 1.5*yystick1:stepy:yystick2*1.5;
      if length(ytickarray)>=L1 && length(ytickarray)<=L2
            ym = ytickarray(round(length(ytickarray)/2));
            ytickarray = ytickarray-ym;
            ytickarray = ytickarray(ytickarray>=yys1 & ytickarray<=yys2 );
            if yys1<0 && ytickarray(1)>=0
                  ytickarray = [ -round((abs(yys1)/2)) ytickarray];
                  ytickarray = unique(ytickarray);
            end
            def = {vect2colon(ytickarray,'Delimiter','off')};
            goags = 0;
      elseif length(ytickarray)>L2
            stepy = stepy*2;
      elseif length(ytickarray)<L1
            stepy = round(stepy/2);
      end
      w=w+1;
end



