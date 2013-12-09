function data = removedc(data, windowsam)

if nargin<2
    windowsam = size(data,2);
end
meanatonset = mean(data(:,1:windowsam),2);
for chr=1:size(data,1)
    data(chr,:) = data(chr,:) - meanatonset(chr);
end
disp('done')