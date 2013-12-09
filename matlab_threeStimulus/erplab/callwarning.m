% callwarning writes the string WARNING in red at the command window
% Examples:
% >>callwarning
%  WARNING: >>
%  
%  >>callwarning('hey')
%  WARNING:hey>>
%  
% >>callwarning('stop ','\n')
%  stop WARNING:
% >>
%
% Author: Javier Lopez-Calderon
function callwarning(varargin)

if nargin<1
        post = '';
        pre  = '';
elseif nargin==1
        post = varargin{1};
        pre = '';
elseif nargin==2
        post = varargin{2};
        pre = varargin{1};
else
        error('callwarning only deals up two inputs')
end
colorw = [0.8 0 0];
try cprintf(colorw, [pre 'WARNING:' post]);catch,fprintf([pre 'WARNING:' post]);end ;


