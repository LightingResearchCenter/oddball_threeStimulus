function plot_debugCOI_onWaveletTransform(time, freq, period, power, coi, handles)

    %% DEBUG
    debugMatFileName = 'tempCoiDebugPlot.mat';
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
    
    %    COI = if specified, then return the Cone-of-Influence, which is a vector
    %        of N points that contains the maximum period of useful information
    %        at that particular time.
    %        Periods greater than this are subject to edge effects.
    
    % Defined as: (analyze_waveletTransform.m)
    % fourier_factor = (4*pi)/(k0 + sqrt(2 + k0^2)); % Scale-->Fourier [Sec.3h
    %       For 'MORLET', k0 (wavenumber), default is 6.        
    % coi = fourier_factor/sqrt(2);
    % coi = coi*dt*[1E-5,1:((n1+1)/2-1),fliplr((1:(n1/2-1))),1E-5];  % COI [Sec.3g]
    %       where dt = 1/fs 
    
    %        This can be used to plot COI lines on a contour plot by doing:
    %
    %              contour(time,log(period),log(power))
    %              plot(time,log(coi),'k')    
    
    scrsz = handles.style.scrsz;
    fig = figure('Color', 'w');
        set(fig, 'Position', [0.4*scrsz(3) 0.05*scrsz(4) 0.45*scrsz(3) 0.91*scrsz(4)])        

        rows = 3;
        cols = 1;
        time = time * 1000; % to ms
    
    contourLevels = 256;
    
    ind = 1;
    sp(ind) = subplot(rows,cols,ind);
    
        hold on
        contourf(time,log(1./freq),log(power), contourLevels, 'EdgeColor', 'none')
        p(ind) = plot(time,log(coi),'k');
        
        lab(ind,1) = xlabel('Time [ms]'); 
        lab(ind,2) = ylabel('log(Period [s])');
        tit(ind) = title(['Wavelet Transform: Single Epoch with COI']);        
        colorbar
        drawnow

    ind = 2;
    sp(ind) = subplot(rows,cols,ind);
    
        hold on
        contourf(time,(freq),log(power), contourLevels, 'EdgeColor', 'none')
        p(ind) = plot(time,(1./coi),'k');

        lab(ind,1) = xlabel('Time [ms]'); 
        lab(ind,2) = ylabel('Frequency [Hz]');                    
        tit(ind) = title(['Wavelet Transform: Single Epoch with COI']);
        colorbar
        ylim([min(freq) max(freq)])
        drawnow
        
    ind = 3;
    sp(ind) = subplot(rows,cols,ind);    
        
        power = analyze_rejectWT_underCOI(time, freq, power, coi, handles);
    
        hold on
        contourf(time,(freq),log(power), contourLevels, 'EdgeColor', 'none')
        p(ind) = plot(time,(1./coi),'k');       
        
        lab(ind,1) = xlabel('Time [ms]'); 
        lab(ind,2) = ylabel('Frequency [Hz]');   
        titStr = sprintf('%s\n%s', ['Wavelet Transform: Single Epoch with COI'], ['Edge effects rejected (under the COI)']);
        tit(ind) = title(titStr);
        ylim([min(freq) max(freq)])
        colorbar
        drawnow
        
    set(p, 'LineWidth', 2)

    set(sp, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-2) 
    set(sp, 'XLim', [min(time) max(time)])
    set(tit, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase, 'FontWeight', 'bold') 
    set(lab, 'FontName', handles.style.fontName, 'FontSize', handles.style.fontSizeBase-1, 'FontWeight', 'bold')
    
    % AUTO-SAVE FIGURE
        try
            if handles.figureOut.debugON == 1 
                drawnow
                dateStr = plot_getDateString(); % get current date as string                
                fileNameOut = ['timeFreqCOI_debug_', dateStr];
                disp(['         ... saving figure to disk (', fileNameOut, '.png]'])
                export_fig(fullfile(handles.path.figuresOut, fileNameOut), handles.figureOut.resolution, handles.figureOut.antialiasLevel, fig)                
            end
        catch err
            err
        end  
        