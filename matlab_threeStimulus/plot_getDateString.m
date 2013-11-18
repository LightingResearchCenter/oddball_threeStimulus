%% DOCUMENTATION
% All the lines starting with %DOC% in similar sense as Java's javadoc

    %DOC%<desc>Main function to do all the tasks needed to import session files<desc/>

    %DOC%<h2>Description</h2>
    %DOC%handles = main_import(handles), call with the existing handles structure and then the function will return you the modified handles with the imported data and settings                

%% CODE
function dateStr = plot_getDateString()

    cl = fix(clock);    
    
    year = num2str(cl(1));
    month = cl(2);
    day = cl(3);
    
    % fix months
    if month < 10; 
        month = ['0', num2str(month)]; 
    else
        month = num2str(month);
    end
    
    % fix days
    if day < 10; 
        day = ['0', num2str(day)]; 
    else
        day = num2str(day);
    end
    
    % return the handle to the text
    dateStr = sprintf('%s%s%s', year, month, day);