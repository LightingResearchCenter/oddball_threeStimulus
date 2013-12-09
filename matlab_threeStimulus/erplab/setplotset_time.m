plotset.ptime.binArray   = binArray;
plotset.ptime.chanArray  = chanArray;
plotset.ptime.chanArray_MGFP = chanArray_MGFP;
plotset.ptime.blcorr     = blcorr;
plotset.ptime.xscale     = xxscale;
plotset.ptime.yscale     = yyscale;
plotset.ptime.linewidth  = linewidth;
plotset.ptime.isiy       = isiy;
plotset.ptime.fschan     = fschan;
plotset.ptime.fslege     = fslege;
%plotset.ptime.meap       = meap;
plotset.ptime.pstyle     = pstyle;
plotset.ptime.errorstd   = errorstd;
plotset.ptime.pbox       = pbox;
plotset.ptime.counterwin = counterwin;
plotset.ptime.holdch     = holdch;
plotset.ptime.yauto      = yauto;
plotset.ptime.yautoticks = yautoticks;
plotset.ptime.xautoticks = xautoticks;
plotset.ptime.binleg     = binleg;
plotset.ptime.chanleg    = chanleg;
plotset.ptime.isMGFP     = isMGFP;
plotset.ptime.ismaxim    = ismaxim;
%plotset.ptime.istopo     = istopo;
plotset.ptime.axsize     = axsize;
plotset.ptime.minorticks = minorticks;
plotset.ptime.linespec   = linespeci;
plotset.ptime.legepos    = legepos;
plotset.ptime.posgui     = posgui;
plotset.ptime.posfig     = posfig;

% adjusting
plotset.ptime.binArray  = plotset.ptime.binArray(plotset.ptime.binArray<=ERP.nbin);
plotset.ptime.chanArray = plotset.ptime.chanArray(plotset.ptime.chanArray<=ERP.nchan);
plotset.ptime.chanArray_MGFP = plotset.ptime.chanArray_MGFP(plotset.ptime.chanArray_MGFP<=ERP.nchan);

% if plotset.ptime.xscale(1) < ERP.xmin*1000
%       plotset.ptime.xscale(1) = ceil(ERP.xmin*1000);
% end
% if plotset.ptime.xscale(2) > ERP.xmax*1000
%       plotset.ptime.xscale(2) = floor(ERP.xmax*1000);
% end