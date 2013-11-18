function p = plot_errorShade(x,mean, std, alpha, acolor)

% Modified from
% http://www.mathworks.com/matlabcentral/fileexchange/29534-stdshade
    % smusall 2010/4/23
    
% and from: 
% http://www.mathworks.com/matlabcentral/fileexchange/26311-shadederrorbar

    % - acolor defines the used color (default is red) 
    % - alpha defines transparency of the shading (default is no shading and black mean line)

% get the limits
uE = mean'+std';
lE = mean'-std'; 

%Make the cordinats for the fill patch
yP=[lE,fliplr(uE)];
xP=[x',fliplr(x')];
    
if exist('alpha','var')==0 || isempty(alpha) 
    p(1:2) = fill(xP, yP, acolor, 'linestyle','none');
    acolor='k';
else
    p(1:2) = fill(xP, yP, acolor, 'FaceAlpha', alpha, 'linestyle', 'none');
end
set(gcf,'renderer','openGL')

if ishold==0
    check=true; 
else
    check=false;
end

hold on;
p(3) = plot(x, mean, 'Color', acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end



