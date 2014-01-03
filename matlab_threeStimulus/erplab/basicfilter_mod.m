% - This function is part of ERPLAB Toolbox -
%
%     EEG         - input dataset
%     chanArray   - channel(s) to filter
%     locutoff    - lower edge of the frequency pass band (Hz)  {0 -> lowpass}
%     hicutoff    - higher edge of the frequency pass band (Hz) {0 -> highpass}
%     filterorder - length of the filter in points {default 3*fix(srate/locutoff)}
%     typef       - type of filter: 0=means IIR Butterworth;  1 = means FIR
%
%     Outputs:
%     EEG         - output dataset
%
% Author: Javier Lopez-Calderon & Steven Luck
% Center for Mind and Brain
% University of California, Davis,
% Davis, CA
% 2009

%b8d3721ed219e65100184c6b95db209bb8d3721ed219e65100184c6b95db209b
%
% ERPLAB Toolbox
% Copyright � 2007 The Regents of the University of California
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

% Modified Petteri Teikari, removed the boundaries and dc_remove option

function [dataOut, ferror, b, a, v, frec3dB_final, xdB_at_fx, orderx] = basicfilter_mod(dataIn, chanArray, Fs, locutoff, hicutoff, filterorder, typef)

    % dataOut = dataIn;
    
        % b  = [bl; bh];
        % a  = [al; ah];

    ferror = 0; % no filter error by default

    if nargin < 1
          help basicfilter
          return
    end
    if exist('filtfilt','file') ~= 2
          error('ERPLAB says: error at basicfilter(). Cannot find the Signal Processing Toolbox');
    end
    if isempty(dataIn)
          error('ERPLAB says: error at basicfilter(). Cannot filter an empty dataset')
    end
    if nargin < 7
          error('ERPLAB says: error at basicfilter(). Please, enter all arguments!')
    end
    if length(dataIn) <= 3*filterorder
        length(dataIn)
        3*filterorder
          msgboxText{1} =  'The length of the data must be more than three times the filter order.';
          title = 'ERPLAB: basicfilter(), filtfilt constraint';
          errorfound(msgboxText, title);
          return
    end
    if locutoff == 0 && hicutoff == 0,
          error('ERPLAB says: Hey dude, I need both lower and higher cuttof in order to filter your data...');
    end

    chanArray = unique(chanArray); % does not allow repeated channels
    fnyquist  = 0.5*Fs;       % half sample rate
    pnts      = length(dataIn);
    numchan   = length(chanArray);
    
    if locutoff >= fnyquist
          error('ERPLAB says: error at basicfilter(). Low cutoff frequency cannot be >= srate/2');
    end
    if hicutoff >= fnyquist
          error('ERPLAB says: error at basicfilter(). High cutoff frequency cannot be >= srate/2');

    end
    if typef>0 && filterorder*3 > pnts  % filtfilt restriction
          fprintf('basicfilter: filter order too high');
          3*filterorder
          length(dataIn)
          error('ERPLAB says: error at basicfilter(). Number of samples must be, at least, 3 times the filter order.');
    end

    % Call the subfunction to "design the filter"
    [b, a, labelf, v, frec3dB_final, xdB_at_fx, orderx] = filter_tf(typef, filterorder, hicutoff, locutoff, Fs);

    if ~v  % something is wrong or turned off
          msgboxText{1} =  'Wrong parameters for filtering.';
          title = 'ERPLAB: basicfilter() error';
          errorfound(msgboxText, title);
          return
    end
    
    % Warning off    
    warning off MATLAB:singularMatrix

    % Filter the whole input data rather than segment-by-segment
    % (trial-by-trial) as done by ERPLAB
    % fprintf('%s filtering input data (fc = %s Hz), please wait...\n\n', labelf, vect2colon(nonzeros([locutoff hicutoff])))
    
    if length(dataIn) > 3*filterorder
        
        if size(b,1)>1            
            if strcmpi(labelf,'Band-Pass')

                % Butterworth bandpass (cascade)                
                dataOut(:,chanArray) = filtfilt(b(1,:),a(1,:), dataIn(:,chanArray));                
                dataOut(:,chanArray) = filtfilt(b(2,:),a(2,:), dataOut(:,chanArray));
          
            else
                
                % Butterworth Notch (parallel)                
                datalowpass   = filtfilt(b(1,:),a(1,:), dataIn(:,chanArray));
                datahighpass  = filtfilt(b(2,:),a(2,:), datalowpass(:,chanArray));                
                dataOut(:,chanArray) = datalowpass + datahighpass;
            end
        else
            % Butterworth lowpass)
            % Butterworth highpass
            % FIR lowpass
            % FIR highpass
            % FIR bandpass
            % FIR notch
            % Parks-McClellan Notch
            dataOut(:,chanArray) = filtfilt(b,a, dataIn(:,chanArray));            
        end
    else
        fprintf('WARNING: EEG segment from sample %d to %d was not filtered.\n', bp1,bp2);
        fprintf('More than 3*filterorder points are required, at least.\n\n');
    end

    fproblems = nnz(isnan(dataIn(:,chanArray)));

    if fproblems > 0
        ferror = 1;    
    end
    %{
          msgboxText = ['Oops! filter is not working properly.\n Data have undefined numerical results.\n'...
                       'We strongly recommend that you change some filter parameters, for instance, decrease filter order.'];
          title = 'ERPLAB: basicfilter() error: undefined numerical results';
          error(sprintf(msgboxText), title);
          return
    end            
    %}

    
    % fprintf('\n')



