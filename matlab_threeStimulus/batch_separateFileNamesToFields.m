function fieldsStructure = batch_separateFileNamesToFields(fileNames)
    
    for i = 1 : length(fileNames)
        
        fileName = fileNames{i};
        
        tmpScan = textscan(fileName, '%s %d %d %s', 'delimiter', '_');

        fieldsStructure{i}.fileName = fileName;
        fieldsStructure{i}.subject = cell2mat(tmpScan{1});
        fieldsStructure{i}.session = tmpScan{2};
        fieldsStructure{i}.intensity = tmpScan{3};
        fieldsStructure{i}.matType = strrep(tmpScan{4}, '.mat', '');
        

        %{
        a = fieldsStructure{i}.fileName
        b = fieldsStructure{i}.subject
        c = fieldsStructure{i}.session
        d = fieldsStructure{i}.intensity
        %}       

    end