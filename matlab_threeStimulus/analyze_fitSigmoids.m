function paramOut = analyze_fitSigmoids(x, y, i, j, k, m, parameters, handles)

    debugMatFileName = 'tempSigmoidFit.mat';
    if nargin == 0
        load('debugPath.mat')
        load(fullfile(path.debugMATs, debugMatFileName))
        close all
    else
        if handles.flags.saveDebugMATs == 1
            path = handles.path;
            save('debugPath.mat', 'path')
            save(fullfile(path.debugMATs, debugMatFileName))            
        end
    end
    
    %x
    %y
    %whos
    
    options = statset('Robust','off');    
    init0 = sigmoid_initCoeffs(x,y); % init values 
    try
        paramOut = nlinfit(x, y, 'sigmoid_4param', init0, options); 
    catch err
        warning(['"', err.identifier, '", not good enough input data, sigmoid not fitted, NaN parameters returned'])
        paramOut = [NaN NaN NaN NaN];
    end
    
    
    function init_params = sigmoid_initCoeffs(x,y)

        % INIT_COEFFS Function to generate the initial parameters for the 4
        % parameter dose response curve.

        parms    = ones(1,4);
        parms(1) = min(y);
        parms(2) = max(y);
        parms(3) = (min(x)+max(x))/2;

        % get input sizes
        sizey    = size(y);
        sizex    = size(x);

        % further fixing
        %{
        if (y(1)-y(sizey)) ./ (x(2)-x(sizex))>0
            parms(4)=(y(1)-y(sizey))./(x(2)-x(sizex));
        else
            parms(4)=1;
        end    
        %}
        parms(4) = 1;
        init_params=parms;
    
    
   