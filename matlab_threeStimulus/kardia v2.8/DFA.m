%% ---- DFA
function [alpha,n,Fn]=DFA(y,varargin)

    % set default values for input arguments
    sliding=0;
    graph=0;
    minbox=4;
    maxbox=floor(length(y)/4);

    % check input arguments
    nbIn = nargin;
    if nbIn > 1
        if ~ischar(varargin{1})
            minbox = varargin{1};
            if ~ischar(varargin{2})
                maxbox = varargin{2};
            else
                error('Input argument missing.');
            end
        end
        for i=1:nbIn-1
            if isequal (varargin{i},'plot'), graph='plot';end
            if isequal (varargin{i},'s'), sliding='s';end
        end
    end

    if nbIn > 5
        error('Too many input arguments.');
    end

    % initialize output variables
    alpha=[];
    n=NaN(1,maxbox-minbox+1);
    Fn=NaN(1,maxbox-minbox+1);

    % transpose data vector if necessary
    s=size(y);
    if s(1)>1
        y=y';
    end

    % substract mean
    y=y-mean(y);

    % integrate time series
    y=cumsum(y);

    N=length(y); % length of data vector

    % error message when box size exceeds permited limits
    if minbox<4 || maxbox>N/4
        disp([mfilename ': either minbox too small or maxbox too large!']);
        return
    end

    % begin loop to change box size
    count=1;
    for n=minbox:maxbox;
        i=1;
        r=N;
        m=[];
        l=[];

        % begin loop to create a new detrended time series using boxes of size n starting
        % from the beginning of the original time series
        while i+n-1<=N % create box size n
            x=y(i:i+n-1);
            x=detrend(x); % linear detrending
            m=[m x];
            if strcmp(sliding,'s')
                i=i+1; % sliding window
            else
                i=i+n; % non-overlapping windows
            end
        end

        % begin loop to create a new detrended time series with boxes of size n starting
        % from the end of the original time series
        while r-n+1>=1
            z=y(r:-1:r-n+1);
            z=detrend(z);
            l=[l z];
            if strcmp(sliding,'s')
                r=r-1;
            else
                r=r-n;
            end
        end

        % calculate the root-mean-square fluctuation of the new time series
        k=[m l]; % concatenate the two detrended time series
        k=k.^2;
        k=mean(k);
        k=sqrt(k);
        Fn(count)=k;
        count=count+1;
    end

    n=minbox:maxbox;

    % plot the DFA
    if strcmp (graph,'plot');
        figure;
        plot(log10(n),log10(Fn))
        xlabel('log(n)')
        ylabel('log(Fn)')
        title('Detrended Fluctuation Analysis')
    end

    % calculate scaling factor alpha
    coeffs    = polyfit(log10(n),log10(Fn),1);
    alpha     = coeffs(1);

end