%
% Author: Javier Lopez-Calderon
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009

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
function button = errorfound(message, title, bkgrncolor, showfig)
if nargin<4
        showfig = 1; % 1 = yes, show fig on error message
end
if nargin<3
        bkgrncolor = [1 0 0];
end
if nargin<2
        title = '???';
end
if nargin<1
        message = 'Wassup?';
end
p = which('eegplugin_erplab');
p = p(1:strfind(p,'eegplugin_erplab.m')-1);
figdir  = fullfile(p,'images');

%
% Get image indices
%
figdir1 = dir(figdir);
ferrorn = regexp({figdir1.name}, 'logoerplaberror\d*.jpg','match', 'ignorecase');
fnames  = cellstr(char([ferrorn{:}]));  
%for i=1:7; fnames{i} = sprintf('logoerplaberror%g.jpg',i);end
try
        [IconData IconCMap]= loadrandimage(fnames);
catch
        [IconData IconCMap]= imread('logoerplaberror1.jpg');
end

button = errorGUI(message, title, IconData, IconCMap, bkgrncolor, showfig);
return


