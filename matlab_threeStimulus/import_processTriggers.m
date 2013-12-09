function triggers = import_processTriggers(triggersRaw, triggerSignals, triggerPrecision, handles)

    debugMatFileName = 'tempTriggerProcessing.mat';
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
      
    % see, http://www.mathworks.com/help/comm/ref/de2bi.html, from Communications System Toolbox
    % download: http://seed.ucsd.edu/mediawiki/images/6/66/De2bi.m          
    
    triggerBinVector = logical(import_De2bi(triggersRaw, triggerPrecision));
    %toc  % -> 1.502745 seconds.
    
        % Alternative ways
        % http://stackoverflow.com/questions/5744576/convert-decimal-to-binary-vector
        
        %{
        tic
        out = zeros(length(triggersRaw), triggerPrecision);
        for rows = 1 : length(triggersRaw)
            out(rows,:) = binary2vector(triggersRaw(rows), triggerPrecision);
        end    
        toc  % -> 69.399970 seconds.

        tic
        out = zeros(length(triggersRaw), triggerPrecision);
        for rows = 1 : length(triggersRaw)
            binVec(rows,:) = dec2bin(triggersRaw(rows), triggerPrecision)-'0';
        end        
        toc  % -> long?

        tic
        binVec = zeros(length(triggersRaw), triggerPrecision);
        for rows = 1 : length(triggersRaw)
            binVec(rows,:) = bitget(triggersRaw(rows),1:triggerPrecision);
        end  
        toc % -> long?
        %}
    
    if handles.flags.showDebugMessages == 1
        [rowsIn, colsIn] = size(triggerBinVector);
        for jjjj = 1 : colsIn    
            uniq = unique(triggerBinVector(:,jjjj));
            if length(uniq) > 1
               disp(['      .. Found trigger: col (channel) index = ', num2str(jjjj)])
            end
        end
    end
           
        
    %% create a structure based on the triggers wanted
            
        % combine both button channels as it doesn't matter which button was
        % used and individual subject should have pressed only either of them
        % LOW is now BUTTON pressed
        triggers.button = logical( ~(triggerBinVector(:,triggerSignals.buttons(1))) + (~triggerBinVector(:,triggerSignals.buttons(2))) );

        triggers.stdTone = triggerBinVector(:,triggerSignals.stdTone);
        triggers.oddTone = triggerBinVector(:,triggerSignals.oddTone);
        triggers.distracter = triggerBinVector(:,triggerSignals.distracterTone);
        
        triggers.audio = triggerBinVector(:,triggerSignals.audioON);
            triggers.audio = import_correctAudioTrigger(triggers.audio, handles);
        
        triggers.recON = triggerBinVector(:,triggerSignals.recON);       
                
      
    %% DEBUG PLOT       
    
        if handles.flags.showDebugPlots == 1
            
            close all
            scrsz = get(0,'ScreenSize'); % get screen size for plotting
            fig = figure('Color', 'w');
                set(fig, 'Position', [0.5*scrsz(3) 0.05*scrsz(4) 0.42*scrsz(3) 0.92*scrsz(4)])
                rows = 5; cols = 1;

            fS = handles.parameters.EEG.srate;
            t = linspace(1,length(triggers.button),length(triggers.button)) / fS;           

            subplot(rows,cols,1)
                area(t, triggers.button); title('Button');
            subplot(rows,cols,2)
                hold on 
                area(t, triggers.oddTone, 'FaceColor', 'r', 'EdgeColor', 'none'); 
                area(t, triggers.stdTone, 'FaceColor', 'k', 'EdgeColor', 'none');
                title('Odd/Std Tone'); legend('Odd', 'Std')            
                hold off
            subplot(rows,cols,3)
                area(t, triggers.irrCycle); title('Irregular cycle');
            subplot(rows,cols,4)
                area(t, triggers.recON); title('REC ON');
            subplot(rows,cols,5)
                area(t, triggersRaw); title('RAW');
                ylim([min(triggersRaw)*0.99 max(triggersRaw)*1.01])
        end
        
        
    % From: http://www.biosemi.com/faq/trigger_signals.htm
    
        %{
        Bit 00 (LSB)
        Trigger Input 1 (High = trigger on)
        Bit 01
        Trigger Input 2 (High = trigger on)
        Bit 02
        Trigger Input 3 (High = trigger on)
        Bit 03
        Trigger Input 4 (High = trigger on)
        Bit 04
        Trigger Input 5 (High = trigger on)
        Bit 05
        Trigger Input 6 (High = trigger on)
        Bit 06
        Trigger Input 7 (High = trigger on)
        Bit 07
        Trigger Input 8 (High = trigger on)
        Bit 08
        Trigger Input 9 (High = trigger on)
        Bit 09
        Trigger Input 10 (High = trigger on)
        Bit 10
        Trigger Input 11 (High = trigger on)
        Bit 11
        Trigger Input 12 (High = trigger on)
        Bit 12
        Trigger Input 13 (High = trigger on)
        Bit 13
        Trigger Input 14 (High = trigger on)
        Bit 14
        Trigger Input 15 (High = trigger on)
        Bit 15
        Trigger Input 16 (High = trigger on)
        Bit 16	
        High when new Epoch is started
        Bit 17	
        Speed bit 0
        Bit 18	
        Speed bit 1
        Bit 19	
        Speed bit 2
        Bit 20	
        High when CMS is within range
        Bit 21	
        Speed bit 3
        Bit 22	
        High when battery is low
        Bit 23 (MSB)	
        High if ActiveTwo MK2
        %}    
    
    
    function out = binary2vector(data,nBits)

        powOf2 = 2.^[0:nBits-1];

        %# do a tiny bit of error-checking
        if data > sum(powOf2)
           error('not enough bits to represent the data')
        end

        out = false(1,nBits);

        ct = nBits;

        while data>0
        if data >= powOf2(ct)
        data = data-powOf2(ct);
        out(ct) = true;
        end
        ct = ct - 1;
        end