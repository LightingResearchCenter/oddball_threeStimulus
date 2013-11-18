function y = pre_nyquistInterpolation(x,ny,dim,nw)
%INTERPFTW 1-D interpolation using FFT method.
%
%   Identical to INTERPFT, but takes optional fourth parameter to 
%   specify aliasing window.
%
%   Y = INTERPFT(X,N,DIM,W) resamples X (along dimension DIM) to 
%   produce interpolated vector Y of length N. Interpolated vector
%   will come from Nyquist window W.
%
%   If W = 0, then INTERPFTW is identical to INTERPFT.
%
%   Class support for data input x:
%      float: double, single
%  
%   See also INTERP1, INTERPFTW.

%   Robert Piche, Tampere University of Technology, 10/93.
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 5.15.4.3 $  $Date: 2005/06/21 19:35:37 $
%
%   Modified by Theodore P. Pavlic, 06/2008.
%   Added Nyquist aliasing support.

error(nargchk(2,4,nargin));

if nargin==2,
  [x,nshifts] = shiftdim(x);
  if isscalar(x), nshifts = 1; end % Return a row for a scalar
  nw = 0;
elseif nargin>3,
  perm = [dim:max(length(size(x)),dim) 1:dim-1];
  x = permute(x,perm);
  if nargin==3, nw = 0; end
end

siz = size(x);
[m,n] = size(x);
if ~isscalar(ny) 
  error('MATLAB:interpft:NonScalarN', 'N must be a scalar.'); 
end

%  If necessary, increase ny by an integer multiple to make ny > m.
if ny > m
   incr = 1;
else
   if ny==0, y=[]; return, end
   incr = floor(m/ny) + 1;
   ny = incr*ny;
end
a = fft(x,[],1);
nyqst = ceil((m+1)/2);
if rem(nw,2) == 1
    a = fftshift(a);
end
b = [zeros(nw*nyqst,n) ; ...
     a(1:nyqst,:) ; ...
     zeros(ny-m-2*nw*nyqst,n) ; ...
     a(nyqst+1:m,:) ; ...
     zeros(nw*nyqst,n)];
if rem(m,2) == 0
   b((nw+1)*nyqst,:) = b((nw+1)*nyqst,:)/2;
   b((1-nw)*nyqst+ny-m,:) = b((nw+1)*nyqst,:);
end
y = ifft(b,[],1);
if isreal(x), y = real(y); end
y = y * ny / m;
y = y(1:incr:ny,:);  % Skip over extra points when oldny <= m.

if nargin==2,
  y = reshape(y,[ones(1,nshifts) size(y,1) siz(2:end)]);
elseif nargin>2,
  y = ipermute(y,perm);
end
