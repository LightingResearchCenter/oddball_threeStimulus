function [importedAlready, matFileToLoad] = import_checkIfAlreadyImported(fileNameIn, matFilesInputFolder)

    matFileToLoad = strrep(fileNameIn, '.bdf', '.mat');

    dirOutputAll = dir(fullfile(matFilesInputFolder, '*.mat')); % Specifies the type of files to be listed        
    dirOutputSpecificFile = dir(fullfile(matFilesInputFolder,matFileToLoad)); % Specifies the type of files to be listed  
    
    importedAlready = ~isempty(dirOutputSpecificFile);
    
    
    