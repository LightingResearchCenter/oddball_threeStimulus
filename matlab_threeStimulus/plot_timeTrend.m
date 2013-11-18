function plot_timeTrend(ERP_Jongsma, epochs_Jongsma, epochs_Jongsma_CNV, epochs_Jongsma_filt, epochs_Jongsma_CNV_filt, ...
                                    epochs_raw, epochs_filt, epochs_CNV_filt, epochs_ep, epochs_ep_CNV, ...
                                    chName, fieldName, cycle, chForIndexERP, chForIndexCNV, ...
                                    contourMode, erpOrCNV, parameters, style, oddballTask, handles)
                                
    debugMatFileName = 'tempTimeTrend.mat';
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
    
    here = 1