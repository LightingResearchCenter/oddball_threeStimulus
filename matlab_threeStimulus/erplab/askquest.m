% - This function is part of ERPLAB Toolbox -
%

function button = askquest(question, tittle, defresp, buttonA, buttonB)

button = '';
if nargin<5
      buttonB = 'NO';
end
if nargin<4
      buttonA = 'Yes';
end
if nargin<3
        defresp = buttonA;
end
%disp(char(question))

TXCOLOR = [1 0.9 0.3];
oldcolor = get(0,'DefaultUicontrolBackgroundColor');
set(0,'DefaultUicontrolBackgroundColor',TXCOLOR)
button = questdlg(question, tittle, buttonA, buttonB, defresp);
set(0,'DefaultUicontrolBackgroundColor',oldcolor)

