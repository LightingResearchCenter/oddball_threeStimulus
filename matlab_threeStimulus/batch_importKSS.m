function KSS = batch_importKSS(inputPath, handles)    
    
    dirOutput     = dir(fullfile(inputPath, 'kss*.txt')); % Specifies the type of files to be listed        
    fileNames     = {dirOutput.name}'; % Prints the filenames from the input folder into a variable            
    
    numberOfFilesFound = length(fileNames);
    
    for file = 1 : numberOfFilesFound
       
        % Get fields in
        fileName = fileNames{file};
        tmpScan = textscan(fileName, '%s %s %d', 'delimiter', '_');       
        
        fileNamesOut{file} = fileName;
        subject{file} = cell2mat(tmpScan{2});
        intensityField{file} = tmpScan{3};
        
        % check intensity
        if intensityField{file} == 0
            intensity = 'dark';
        elseif intensityField{file} == 10
            intensity = 'dim';
        elseif intensityField{file} == 40
            intensity = 'bright';
        else            
            a = intensityField{file}
            error('Error in the intensity field!')
        end
        
        importTemp = importdata(fullfile(inputPath, fileNamesOut{file}));
        KSS.(intensity).(subject{file}) = importTemp;
        
    end
    