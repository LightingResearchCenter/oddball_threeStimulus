function [hMRAplotAxes, hXplotAxes] = plot_imodwt_mra_mod(DJt, SJ0t, X, mra_att, ...
                                                  title_str, xaxis, xlabel_str, ...
                                                  MRAplotAxesProp, XplotAxesProp, ...
                                                  J0, level_range, plotOpts, ...
                                                  masterPlotFrame, xaxis_range)
% plot_imodwt_mra -- Plot the inverse MODWT multiresolution analysis detail and smooth coefficients and original time series.
%
%****f* wmtsa.dwt/plot_modwt_mra
%
% NAME
%   plot_imodwt_mra -- Plot the inverse MODWT multiresolution analysis detail and smooth coefficients and original time series.
%
% SYNOPSIS
%   [hDplotAxes, hXplotAxes] = plot_imodwt_mra(DJt, SJ0t, [X], mra_att, ... 
%                                   [title_str], [xaxis], [xlabel_str], ...
%                                   [MRAplotAxesProp], [XplotAxesProp], ...
%                                   [J0], [level_range], [plotOpts], ...
%                                   [masterPlotFrame], [xaxis_range])
%   [hDplotAxes, hXplotAxes] = plot_imodwt_mra(DJt,   [], [X], mra_att, ... 
%                                   [title_str], [xaxis], [xlabel_str], ...
%                                   [MRAplotAxesProp], [XplotAxesProp], ...
%                                   [J0], [level_range], [plotOpts], ...
%                                   [masterPlotFrame], [xaxis_range])
%   [hDplotAxes, hXplotAxes] = plot_imodwt_mra( [], SJ0t, [X], mra_att, ... 
%                                   [title_str], [xaxis], [xlabel_str], ...
%                                   [MRAplotAxesProp], [XplotAxesProp], ...
%                                   [J0], [level_range], [plotOpts], ...
%                                   [masterPlotFrame], [xaxis_range])
%
% INPUTS
%   DJt           =  NxJ array of reconstituted detailed data series.
%                    where N = number of time intervals,
%                          J = number of levels
%   SJO           =  Nx1 vector of reconstituted smoothed data series.
%   X             =  (optional) vector of observations of size N.
%   * mra_att        -- MODWT transform attributes (struct).
%   title_str     =  (optional) character string or cell array of strings containing title of plot.
%   xaxis         =  (optional) vector of values to use for plotting x-axis.
%   xlabel_str    =  (optional) character string or cell array of strings containing label x-axis.
%   WplotAxesProp =  (optional) structure containing axes property values to
%                    override for W subplot.
%   XplotAxesProp =  (optional) structure containing axes property values to
%                    override for X subplot.
%   J0            =  (optional) override value of J0, if J ~= J0 or
%                    if max(level_range) ~= J0.
%   level_range   =  (optional) number or range of numbers indicating subset of
%                    levels (j's) to plot.
%   plotOpts      =  (optional) structure containing plot options.
%   masterPlotFrame = (optional) structure containing coordinates to 
%                     place X and W plots.
%   xaxis_range   =  (optional) range of xaxis values to plot.
%
% OUTPUTS
%   hMRAplotAxes  =  (optional) handle to axes for multiresolution analysis (MRA) subplot.
%   hXplotAxes     =  (optional) handle to axes for original data series (X) subplot.
%
% SIDE EFFECTS
%
%
% DESCRIPTION
%   plot_imodwt_mra plots the inverse MODWT detail and smowth coefficients
%   and optionally the original data series.  Either or both the inverse
%   MODWT detail and smooth coefficients may be plotted.  
%
%   By default, the detail coefficients (DJt) and scaling coefficient at
%   level J0 (SJ0t) are plotted.  A subrange of wavelet coefficient levels
%   may be specified by via the parameter, level_range = [lower:upper]
%   
%
% EXAMPLE
%
%
% NOTES
%
%
% REFERENCES
%
%
% SEE ALSO
%   modwt, modwt_filter, overplot_modwt_cir_shift_coef_bdry,
%   multi_yoffset_plot
%

% AUTHOR
%   Charlie Cornish
%
% CREATION DATE
%   2003/05/01
%
% COPYRIGHT
%
%
% CREDITS
%
%
% REVISION
%   $Revision: 612 $
%
%***

% $Id: plot_imodwt_mra.m 612 2005-10-28 21:42:24Z ccornish $

% Initialize constants and other parameters
ylabel_xoffset = .015;

% Master Plot Frame - contains all subplot
default_masterPlotFrame.left = .13;
default_masterPlotFrame.bottom = .11;
default_masterPlotFrame.width = .775;
default_masterPlotFrame.height = .805;
default_masterPlotFrame.subplot_yspacing = .010;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Check Input Arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

usage_str = ['Usage:  [hMRAplotAxes, hXplotAxes] = ', mfilename, ...
             '(DJt, SJ0t, [X], mra_att' ...
             ', [title_str], [xaxis], [xaxis_label]', ...
             ', [MRAplotAxesProp], [XplotAxesProp]' ...
             ', [J0], [level_range], [plotOpts], ', ...
             ', [masterPlotFrame], [xaxis_range])'];

%%  Check input arguments and set defaults.
error(nargerr(mfilename, nargin, [1:14], nargout, [0:2], 1, usage_str, 'struct'));

if (~exist('mra_att', 'var') || isempty(mra_att) )
  error('Must specify the MODWT_MRA attribute structure.');
end  

wtfname = mra_att.WTF;
NX  = mra_att.NX;
NW = mra_att.NW;
J0 = mra_att.J0;
boundary = mra_att.Boundary;


%% Get a valid wavelet transform filter coefficients struct.
if (ischar(wtfname))
  try
    [wtf_s] = modwt_filter(wtfname);
  catch
    rethrow(lasterror);
  end
else
  error('WMTSA:invalidWaveletTransformFilter', ...
        encode_errmsg('WMTSA:invalidWaveletTransformFilter', wmtsa_err_table, 'wtf'));
end
  
wtfname = wtf_s.Name;
gt = wtf_s.g;
ht = wtf_s.h;

% Initialize and set the plot flags
plot_DJt = 0;
plot_SJ0t = 0;
plot_X = 0;
plot_modwt_boundary = 1;

if (exist('DJt', 'var') && ~isempty(DJt) )
  plot_DJt = 1;
end  

if (exist('SJ0t', 'var') && ~isempty(SJ0t) )
  plot_SJ0t = 1;
end

if (exist('X', 'var') && ~isempty(X) )
  plot_X = 1;
end

if (plot_DJt == 0 && plot_SJ0t == 0)
  error('Must specify either DJt or SJ0t, or both for plotting');
end

if (~exist('plotOpts', 'var'))
  plotOpts = [];
end

if (isfield(plotOpts, 'PlotDJt'))
  if (plotOpts.PlotDJt)
    plot_DJt = 1;
  else
    plot_DJt = 0;
  end
end

if (isfield(plotOpts, 'PlotSJ0t'))
  if (plotOpts.PlotSJ0t)
    plot_SJ0t = 1;
  else
    plot_SJ0t = 0;
  end
end

if (isfield(plotOpts, 'PlotX'))
  if (plotOpts.PlotX)
    plot_X = 1;
  else
    plot_X = 0;
  end
end

if (plot_DJt == 0 && plot_SJ0t == 0)
  error('Must specify either DJt or SJ0t, or both for plotting');
end

if (isfield(plotOpts, 'PlotMODWTBoundary'))
  if (plotOpts.PlotMODWTBoundary)
    plot_modwt_boundary = 1;
  else
    plot_modwt_boundary = 0;
  end    
end

if (~exist('level_range', 'var') || isempty(level_range))
  level_range = [];
else
  error(argterr(mfilename, level_range, 'posint', [], 1, '', 'struct'));
end

% Initialized and override masterPlotFrame

if (~exist('masterPlotFrame', 'var') || isempty(masterPlotFrame))
  masterPlotFrame = default_masterPlotFrame;
else
  if (~isfield(masterPlotFrame, 'left'))
    masterPlotFrame.left = default_masterPlotFrame.left;
  end
  if (~isfield(masterPlotFrame, 'bottom'))
    masterPlotFrame.bottom = default_masterPlotFrame.bottom;
  end
  if (~isfield(masterPlotFrame, 'width'))
    masterPlotFrame.width = default_masterPlotFrame.width;
  end
  if (~isfield(masterPlotFrame, 'height'))
    masterPlotFrame.height = default_masterPlotFrame.height;
  end
  if (~isfield(masterPlotFrame, 'wsubplot_yspacing'))
    masterPlotFrame.subplot_yspacing = default_masterPlotFrame.subplot_yspacing;
  end
end

% Initial default values for the MRA (coefficients) plot frame.
MRAplotFrame.left = masterPlotFrame.left;
MRAplotFrame.bottom = masterPlotFrame.bottom;
MRAplotFrame.width = masterPlotFrame.width;
MRAplotFrame.height = masterPlotFrame.height;
  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Transform data for plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MRA = [];
MRAlabels = {};



jj = 0;

if (plot_DJt)
  if(~exist('level_range', 'var') || isempty(level_range))
    level_range = 1:J0;
  end
  for (j = level_range)
    jj = jj + 1;
    MRA(:,jj) = DJt(:,j);
    Kj=2^j;
    %MRAlabels{jj} = ['D\rm_{', int2str(j), '}'];
    MRAlabels{jj} = ['\Deltat = \rm{', int2str(Kj), '}'];
  end
end

if (plot_SJ0t)
  jj = jj + 1;
  MRA(:,jj) = SJ0t;
  Kj0=2^(J0+1);
  %MRAlabels{jj} = ['S\rm_{', int2str(J0), '}'];
  MRAlabels{jj} = ['\Deltat = +\rm{', int2str(Kj0), '}'];
end

N = size(MRA,1);

%% Truncate for reflection BC's, when NW = 2 * N;
if (size(MRA,1) == NW & NW == 2 * NX)
  N = NX;
  MRA = MRA(1:N,:,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Setup up x-axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If xaxis is not specified, use the sample points
if (~exist('xaxis', 'var') || isempty(xaxis))
  xaxis = (1:N);
end

if (exist('xaxis_range', 'var') && ~isempty(xaxis_range))
  indices = find(xaxis >= xaxis_range(1) & xaxis <= xaxis_range(end));
  xaxis = xaxis(indices);
  MRA = MRA(indices,:);
  X = X(indices,1);
else
  xaxis_range = [];
end

xaxis_xmin = min(xaxis);
xaxis_xmax = max(xaxis);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Parse the Xplot options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plot_X)

  % Check that length(X) == length(MRA)
  if (length(X) ~= length(MRA))
    error(['Length of X must be equal nrows of DJt and SJ0t']);
  end

  % Initialize default values
  Xmin = min(X);
  Xmax = max(X);

  XplotAxesPropName = {};
  XplotAxesPropVal = {};
  
  if (exist('XplotAxesProp', 'var') && ~isempty(XplotAxesProp))

    % Overide default values for Xmin, Xmin (= YLim min and max for Xplot)
    if (isfield(XplotAxesProp, 'YLim'))
      ylim = XplotAxesProp.YLim;
      Xmin = ylim(1);
      Xmax = ylim(2);
    end

    % Populate cell arrays for Xplot axes names and properties
    axes_field_names = fieldnames(XplotAxesProp);
    nfields = length(axes_field_names);
    
    for (i = 1:nfields)
      fname = axes_field_names{i};
      XplotAxesPropName(i) = {fname};
      fval = XplotAxesProp.(fname);
      XplotAxesPropVal(i) = {fval};
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  Override x-axis min, max values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (exist('XplotAxesProp', 'var'))
  % Overide default values for xaxis_xmin, xaxis_xmax from XplotAxesProp
  if (isfield(XplotAxesProp, 'XLim'))
    xlim = XplotAxesProp.XLim;
    xaxis_xmin = xlim(1);
    xaxis_xmax = xlim(2);
  end

elseif (exist('MRAplotAxesProp', 'var'))
  % Overide default values for xaxis_xmin, xaxis_xmax from MRAplotAxesProp
  if (isfield(MRAplotAxesProp, 'XLim'))
    xlim = MRAplotAxesProp.XLim;
    xaxis_xmin = xlim(1);
    xaxis_xmax = xlim(2);
  end
else
  % Do nothing - no override
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Create Xplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (plot_X)

  % In creating the subplots for X and MRA, scale the subplots so that X and MRA
  % are plotted on the same y scale.

  MRAmax = max(MRA);
  MRAmin = min(MRA);
  MRAtotal_range = sum(abs(MRAmax - MRAmin));
  Xtotal_range = abs(Xmax - Xmin);
  Yaxis_total_range = MRAtotal_range + Xtotal_range;
  MRAplot_ypercent = MRAtotal_range / Yaxis_total_range;
  Xplot_ypercent = Xtotal_range / Yaxis_total_range;

  XplotFrame.height = masterPlotFrame.height * Xplot_ypercent;
  MRAplotFrame.height = masterPlotFrame.height * MRAplot_ypercent;

  XplotFrame.left = masterPlotFrame.left;
  XplotFrame.bottom = masterPlotFrame.bottom;
  XplotFrame.width = masterPlotFrame.width;
  
  Position = [XplotFrame.left, ...
              XplotFrame.bottom ...
              XplotFrame.width, ...
              XplotFrame.height];

  % Create a subplot for X;
  axes('Position', Position, ...
       'XAxisLocation', 'bottom', ...
       'YAxisLocation', 'left');

  line(xaxis, X);

  set(gca, XplotAxesPropName, XplotAxesPropVal);
  set(gca, 'XLim', [xaxis_xmin xaxis_xmax]);
  set(gca, 'YLim', [Xmin Xmax]);

  xscale = (xaxis_xmax - xaxis_xmin) / masterPlotFrame.width;
  ylabel_xpos = xaxis_xmax + (xscale * ylabel_xoffset);

  if ( min(X) < 0 && max(X) > 0)
    ylabel_ypos = 0;
  else
    ylabel_ypos = Xmin + (Xmax - Xmin) / 2;
  end
  
  text(ylabel_xpos, ylabel_ypos, '\bf Signal X(t)', ...
       'HorizontalAlignment', 'Left', ...
       'VerticalAlignment', 'Middle');

  if (exist('xlabel_str', 'var') && ~isempty(xlabel_str))
    xlabel(xlabel_str)
  end

  hXplotAxes = gca;
  
  MRAplotFrame.bottom = masterPlotFrame.bottom + XplotFrame.height + masterPlotFrame.subplot_yspacing;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Create MRAplot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MRAplotPosition = [MRAplotFrame.left, ...
                   MRAplotFrame.bottom ...
                   MRAplotFrame.width, ...
                   MRAplotFrame.height];

hMRAplotAxes  = axes('Tag', 'MRAplot',  ...
                     'Position', MRAplotPosition);

% Plot the MRA (details and smooth)

MRAplotAxesProp.XLim = [xaxis_xmin xaxis_xmax];
hMRAplotAxes = multi_yoffset_plot(xaxis, MRA, MRAlabels, MRAplotAxesProp);

% Turn off tick mark labeling for MRA plot.
set(hMRAplotAxes, 'YTickLabel', {});
set(hMRAplotAxes, 'XTickLabel', {});

if (~plot_X)
  % If no Xplot ...
  %   - add xaxis label if string value provided
  if (exist('xlabel_str', 'var') && ~isempty(xlabel_str))
    xlabel(xlabel_str)
  end
  %  - turn on tickmark lablels if overrided.
  if (isfield(MRAplotAxesProp, 'XTickLabel'))
    set(hMRAplotAxes, 'XTickLabel', MRAplotAxesProp.XTickLabel);
  end
  if (isfield(MRAplotAxesProp, 'YTickLabel'))
    disp(MRAplotAxesProp.YTickLabel);
    set(hMRAplotAxes, 'YTickLabel', MRAplotAxesProp.YTickLabel);
  end
end

if (exist('title_str', 'var') && ~isempty(title_str))
  suptitle(title_str)
end

% TODO:  For now, don't plot circular shift boundaries if xaxis_range specified
%        Add at later date
if (plot_modwt_boundary && isempty(xaxis_range))
  overplot_imodwt_cir_mra_bdry(hMRAplotAxes, DJt, SJ0t, mra_att, ...
                                     xaxis, J0, level_range);
end

return
