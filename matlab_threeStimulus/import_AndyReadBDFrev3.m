function [INFO,M,triggers]=import_AndyReadBDFrev3(FileName, SignalsToRead, invertPolarity)
% SignalsToRead is the vector of channels to read, e.g., [1 2 17] reads
% channels 1, 2 and channel 17

    %% HEADER

    bytesPerSample = 3;
    try
        fid=fopen(FileName,'r','ieee-le');          
    catch err
        err
        error(err.identifier)
    end
        
    if fid<0 
        fprintf(2,['Error Loading BDF: File ' FileName ' not found\n']);  
        return;
    end;

    INFO.FileName = FileName;

    H1=char(fread(fid,256,'char')');
    INFO.HeadLen = str2num(H1(185:192));  % 8 Byte  Length of Header
    INFO.DataFormat = char(H1(193:236));  % Version of data format (16 or 24-bit)		
    INFO.NRec = str2num(H1(237:244));     % 8 Byte  # of data records
    INFO.Dur = str2num(H1(245:252));      % 8 Byte  # duration of data record in sec
    INFO.NS = str2num(H1(253:256));       % 8 Byte  # of signals

    INFO.Label = char(fread(fid,[16,INFO.NS],'char')');		
    INFO.Transducer = char(fread(fid,[80,INFO.NS],'char')');	
    INFO.PhysDim = char(fread(fid,[8,INFO.NS],'char')');	

    INFO.PhysMin= str2num(char(fread(fid,[8,INFO.NS],'char')'));
    INFO.PhysMax= str2num(char(fread(fid,[8,INFO.NS],'char')'));
    INFO.DigMin = str2num(char(fread(fid,[8,INFO.NS],'char')'));	
    INFO.DigMax = str2num(char(fread(fid,[8,INFO.NS],'char')'));	
    INFO.PreFilt= char(fread(fid,[80,INFO.NS],'char')');	
    tmp = fread(fid,[8,INFO.NS],'char')';       % Samples per data record
    INFO.SamplesPerRec = str2num(char(tmp));	% Samples per data record

    fseek(fid,32*INFO.NS,0);                    % Position file marker at end of header
            
    INFO.srate = INFO.SamplesPerRec(1); % save the sampling rate
    
    % ***************************** End of Header Section ********************


    if find(diff(INFO.NS)~=0)~=[]
        fprintf(2,['Number of samples for each channel is different \n']);
        M = NaN;
        return
    end

    if max(SignalsToRead) > INFO.NS
        fprintf(2,['Requested number of signals exceeds that in file \n']);
        return
    end

    totalLength = INFO.NRec*max(INFO.SamplesPerRec);
    M = NaN*ones(totalLength,length(SignalsToRead));

    % fseek(fid,(INFO.HeadLen),'bof'); % Should already be at this position
    if INFO.NRec==-1
        fprintf(2,['Number of Records recorded is "unknown" (-1), cannot continue \n']);
    end
    for i=1:INFO.NRec
        for j=1:length(SignalsToRead)
            filePosition = INFO.HeadLen + (i-1)*bytesPerSample*sum(INFO.SamplesPerRec)...
                + bytesPerSample*(sum(INFO.SamplesPerRec(1:SignalsToRead(j)))-INFO.SamplesPerRec(SignalsToRead(j)));
            index = (i-1)*INFO.SamplesPerRec(SignalsToRead(j)) + 1;
            fseek(fid,filePosition,'bof');
            [M(index:index+INFO.SamplesPerRec(SignalsToRead(j))-1,j), count]=fread(fid,INFO.SamplesPerRec(SignalsToRead(j)),'bit24');
        end
    end

    
    
    %% TRIGGER CHANNELS        
    %if SignalsToRead(end)==INFO.NS
        
        triggers = M(:,end);
        % triggers(1:30)
        %triggers = dec2bin(triggers,24);
        %triggersBin(1:30)
        % dec2bin(unique(triggers))
        % plot(triggers)
        %whos
        %pause
        
        % remove the triggers from the EEG recording
        M = M(:,1:end-1);     
        
        % Mask off upper 8 bits of trigger channel if trigger channel is included
        % Ahex = dec2hex(M(:,end),6); % Convert to hexadecimal
        % Ahex = Ahex(:,3:6); % remove upper 8 bits
        % M(:,end) = hex2dec(Ahex); % Convert back to decimal        
        
    % end

    %% Calibration of channels
    
        
        CalFactors = (INFO.PhysMax-INFO.PhysMin)./(INFO.DigMax-INFO.DigMin);
        for i = 1 : length(SignalsToRead) - 1 % not the trigger channel
            
            M(:,i) = M(:,i) .* CalFactors(i);
            if invertPolarity == 1
                M(:,i) = M(:,i) * -1;
            end
            
        end
        
        
        return
