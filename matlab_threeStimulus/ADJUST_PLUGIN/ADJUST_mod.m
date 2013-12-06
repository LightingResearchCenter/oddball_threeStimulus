
% ADJUST() - Automatic EEG artifact Detector 
% with Joint Use of Spatial and Temporal features
%
% Usage:
%   >> [art, horiz, vert, blink, disc,...
%         soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
%         soglia_GDSF, GDSF, soglia_V, nuovaV]=ADJUST(EEG,out);
%
% Inputs:
%   EEG        - current dataset structure or structure array (has to be epoched)
%   out        - (string) report file name 
%
% Outputs:
%   art        - List of artifacted ICs
%   horiz      - List of HEM ICs 
%   vert       - List of VEM ICs   
%   blink      - List of EB ICs     
%   disc       - List of GD ICs     
%   soglia_DV  - SVD threshold      
%   diff_var   - SVD feature values
%   soglia_K   - TK threshold      
%   meanK      - TK feature values
%   soglia_SED - SED threshold      
%   SED        - SED feature values
%   soglia_SAD - SAD threshold      
%   SAD        - SAD feature values
%   soglia_GDSF- GDSF threshold      
%   GDSF       - GDSF feature values
%   soglia_V   - MEV threshold      
%   nuovaV     - MEV feature values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% ADJUST
% Automatic EEG artifact Detector with Joint Use of Spatial and Temporal
% features
% Developed May2007 - October2008
% Andrea Mognon and Marco Buiatti
% CIMeC - Center for Mind/Brain Science, University of Trento
% Last update: 26/11/2009 by Andrea Mognon
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Reference paper:
% Mognon, Jovicich, Bruzzone, Buiatti, ADJUST: An Automatic EEG artifact Detector based on the
% Joint Use of Spatial and Temporal features. Reviewed
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2009 Andrea Mognon and Marco Buiatti, 
% Center for Mind/Brain Sciences, University of Trento, Italy
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% VERSIONS LOG
% 
% V2 (07 OCTOBER 2010) - by Andrea Mognon
% Added input 'nchannels' to compute_SAD and compute_SED_NOnorm;
% this is useful to differentiate the number of ICs (n) and the number of
% sensors (nchannels);
% bug reported by Guido Hesselman on October, 1 2010.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [art, horiz, vert, blink, disc,...
        soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
        soglia_GDSF, GDSF, soglia_V, nuovaV, soglia_D, maxdin]=ADJUST (EEG,out,handles)

    [~, handles.flags] = init_DefaultSettings(); % use a subfunction        
    if handles.flags.saveDebugMATs == 1
        debugMatFileName = 'tempADJUST.mat';
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
    end

    % size(EEG.data)
    %size(EEG.icaact)
    % EEG
    %EEG.icaact = [];
    totalSamplesPerChannel = size(EEG.data, 2) * size(EEG.data, 3);

    %% Settings

        % ----------------------------------------------------
        % |  Change experimental settings in this section    |
        % ----------------------------------------------------

        % ----------------------------------------------------
        % |  Collect useful data from EEG structure          |
        % ----------------------------------------------------

        %number of ICs=size(EEG.icawinv,1);
        %number of time points=size(EEG.data,2);
        if length(size(EEG.data))==3   
            num_epoch=size(EEG.data,3);       
        else    
            num_epoch=0;
        end

        %% Check the presence of ICA activations
        if isempty(EEG.icaact)
            disp('EEG.icaact not present. Recomputed from data.');
            if length(size(EEG.data))==3
        %         EEG.icaact = EEG.icaweights*EEG.icasphere*reshape(EEG.data, size(EEG.icawinv,1), num_epoch*size(EEG.data,2));
        %         EEG.icaact = reshape(EEG.icaact,size(EEG.icawinv,1),size(EEG.data,2), num_epoch);
                 try
                    EEG.icaact = reshape(EEG.icaweights*EEG.icasphere*reshape(EEG.data,[size(EEG.data,1)...
                                         size(EEG.data,2)*size(EEG.data,3)]),[size(EEG.data,1) size(EEG.data,2) size(EEG.data,3)]);
                 catch err
                    err
                    error('the dimension not really function, in the pre_artifactFaster... the unmixing matrix is used to calculate ICA Activations')
                 end
                 
            else
                EEG.icaact = EEG.icaweights*EEG.icasphere*EEG.data;
            end
        end

        topografie=EEG.icawinv'; %computes IC topographies

        %% Topographies and time courses normalization
        disp('          - Normalizing topographies...')
        disp('          - Scaling time courses...')

        for i=1:size(EEG.icawinv,2) % number of ICs    
            ScalingFactor=norm(topografie(i,:));    
            topografie(i,:)=topografie(i,:)/ScalingFactor;

            if length(size(EEG.data))==3
                EEG.icaact(i,:,:)=ScalingFactor*EEG.icaact(i,:,:);
            else
                EEG.icaact(i,:)=ScalingFactor*EEG.icaact(i,:);
            end    
        end

        % Variables memorizing artifacted ICs indexes
        blink=[];
        horiz=[];
        vert=[];
        disc=[];

    %% Feature extraction

        %GDSF - General Discontinuity Spatial Feature
        disp('           - General Discontinuity Spatial Feature...')
        try
            GDSF = compute_GD_feat(topografie, EEG.chanlocs, size(EEG.icawinv,2));
        catch err
            if strcmp(err.identifier, 'MATLAB:badsubscript')
                disp('             .. could not be computed, not enough channels')
                warning(['For educational purposes the script now runs until the end, but when we have no "spatial articats"', ...
                         'no artifacts are found in the end with the spatiotemporal artifact decision'])
                GDSF = zeros(1, size(EEG.icawinv,2));
            else
                err
                error('Some other error here')
            end
        end

        %SED - Spatial Eye Difference
        disp('           - SED - Spatial Eye Difference...')
        try
            [SED,medie_left,medie_right]=computeSED_NOnorm(topografie,EEG.chanlocs,size(EEG.icawinv,2),size(EEG.icawinv,1)); 
        catch err
            if strcmp(err.identifier, 'MATLAB:badsubscript')
                disp('             .. could not be computed, not enough channels')
                SED = zeros(1, size(EEG.icawinv,2));            
                medie_left=zeros(1,size(EEG.icawinv,2));
                medie_right=zeros(1,size(EEG.icawinv,2));
            else
                err
                error('Some other error here')
            end
        end


        %SAD - Spatial Average Difference
        disp('           - SAD - Spatial Average Difference...')
        try
            [SAD,var_front,var_back,mean_front,mean_back]=computeSAD(topografie,EEG.chanlocs,size(EEG.icawinv,2),size(EEG.icawinv,1));
        catch err
            if strcmp(err.identifier, 'MATLAB:badsubscript')
                disp('             .. could not be computed, not enough channels')
                SAD = zeros(1, size(EEG.icawinv,2));            
                var_front=zeros(1,size(EEG.icawinv,2));
                var_back=zeros(1,size(EEG.icawinv,2));
                mean_front=zeros(1,size(EEG.icawinv,2));
                mean_back=zeros(1,size(EEG.icawinv,2));
            else
                err
                error('Some other error here')            
            end
        end


        %SVD - Spatial Variance Difference between front zone and back zone
        diff_var=var_front-var_back;

        disp(' '); disp('             --> Done.'); disp(' ')

        % epoch dynamic range, variance and kurtosis
        K=zeros(num_epoch,size(EEG.icawinv,2)); %kurtosis
        Vmax=zeros(num_epoch,size(EEG.icawinv,2)); %variance
        disp('           - Computing variance and kurtosis of all epochs...')

        if length(size(EEG.icaact)) < 3
            % disp('2D')
            % if epochs are concatenated
            EEG.icaact = pre_epochsVectorToMatrix(EEG.icaact,EEG.pnts);
        else
            disp('3D')
        end
        
        for ic = 1:size(EEG.icawinv,2) % number of ICs
            for ep = 1:num_epoch
                Vmax(ep,ic)=var(EEG.icaact(ic,:,ep));
                K(ep,ic)=kurtosis(EEG.icaact(ic,:,ep));
            end
        end

        % TK - Temporal Kurtosis
        disp('           - Temporal Kurtosis...')
        meanK=zeros(1,size(EEG.icawinv,2));

        for i=1:size(EEG.icawinv,2)
            if num_epoch>100
                meanK(1,i)=trim_and_mean(K(:,i)); 
            else
                meanK(1,i)=mean(K(:,i));
            end

        end

        %MEV - Maximum Epoch Variance
        disp('           - Maximum epoch variance...')

        maxvar=zeros(1,size(EEG.icawinv,2));
        meanvar=zeros(1,size(EEG.icawinv,2));

        for i=1:size(EEG.icawinv,2)
            if num_epoch>100
             maxvar(1,i)=trim_and_max(Vmax(:,i)');
             meanvar(1,i)=trim_and_mean(Vmax(:,i)');
            else 
             maxvar(1,i)=max(Vmax(:,i));
             meanvar(1,i)=mean(Vmax(:,i));
            end
        end

        % MEV in reviewed formulation:
        nuovaV=maxvar./meanvar;



    %% Thresholds computation

        disp('           - Computing EM thresholds...')
        
            % K
            if sum(meanK == 0) == length(meanK)
                disp('              .. meanK could not be computed above, returning NaN')
                soglia_K = NaN;
            else
                soglia_K=EM(meanK);
            end

            % SED
            if sum(SED == 0) == length(SED)
                disp('              .. SED could not be computed above, returning NaN')
                soglia_SED = NaN;
            else
                soglia_SED=EM(SED);
            end

            % SAD
            if sum(SAD == 0) == length(SAD)
                disp('              .. SED could not be computed above, returning NaN')
                soglia_SAD = NaN;
            else
                soglia_SAD=EM(SAD);
            end

            % GDSF
            if sum(GDSF == 0) == length(GDSF)
                disp('              .. GDSF could not be computed above, returning NaN')
                soglia_GDSF = NaN;
            else
                soglia_GDSF=EM(GDSF);    
            end

            % New V
            if sum(nuovaV == 0) == length(nuovaV)
                disp('              .. V could not be computed above, returning NaN')
                soglia_V = NaN;
            else
                soglia_V=EM(nuovaV); 
            end
            

        disp(' '); disp('             --> Done.'); disp(' ')


    %% Output file header

        % ----------------------------------------------------
        % |  Opens report file and writes header             |
        % ----------------------------------------------------

        file=fopen(out,'w');
        fprintf(file,'ADJUST\n');
        fprintf(file,'Automatic EEG artifacts Detector with Joint Use of Spatial and Temporal features\n\n');
        fprintf(file,['Analyzed dataset: ' EEG.filename '\n']);
        fprintf(file,'Analysis carried out on the %d Independent Components\n\n',size(EEG.icawinv,2));

    %% Horizontal eye movements (HEM)

        disp('           - Evaluating Horizontal movements...')

        % ----------------------------------------------------
        % |  Writes HEM header in the report file            |
        % ----------------------------------------------------

        fprintf(file,'> HEM - Horizontal movements\n\n');
        fprintf(file,'Classification based on features:\n\n');
        fprintf(file,'SED - Spatial eye difference (threshold=%f)\n',soglia_SED);
        fprintf(file,'MEV - Maximum epoch variance (threshold=%f)\n\n',soglia_V);
        fprintf(file,'ICs with Horizontal eye movements:\n');

        horiz=intersect(intersect(find(SED>=soglia_SED),find(medie_left.*medie_right<0)),...
            (find(nuovaV>=soglia_V)));

        hor_bool=1; %true if there are artifacted ICs

        if isempty(horiz) %no IC found

            fprintf(file,'/ \n');
            hor_bool=0;

        else

            fprintf(file,[num2str(horiz) '\n']);
            fprintf(file,'\n');

        end


    %% Vertical eye movements (VEM)

        disp('            - Evaluating Vertical movements...')

        % ----------------------------------------------------
        % |  Writes VEM header in the report file            |
        % ----------------------------------------------------

        fprintf(file,'>> VEM - Vertical movements\n\n');
        fprintf(file,'Classification based on features:\n\n');
        fprintf(file,'SAD - Spatial average difference (threshold=%f)\n',soglia_SAD);
        fprintf(file,'MEV - Maximum epoch variance (threshold=%f)\n',soglia_V);
        fprintf(file,'ICs with Vertical eye movements:\n');

        vert=intersect(intersect(find(SAD>=soglia_SAD),find(medie_left.*medie_right>0)),...
            intersect(find(diff_var>0),find(nuovaV>=soglia_V)));

        ver_bool=1; %true if there are artifacted ICs
        if isempty(vert) %no artifact found
            fprintf(file,'/ \n');
            ver_bool=0;
        else
            fprintf(file,[num2str(vert) '\n']);
            fprintf(file,'\n');    
        end


    %% Eye Blink (EB)

        disp('             - Evaluating Blinks...')
        % ----------------------------------------------------
        % |  Writes EB header in the report file             |
        % ----------------------------------------------------
        fprintf(file,'>>> EB - Blinks\n\n');
        fprintf(file,'Classification based on features:\n\n');
        fprintf(file,'SAD (threshold=%f)\n',soglia_SAD);
        fprintf(file,'TK - Temporal kurtosis (threshold=%f)\n',soglia_K);
        fprintf(file,'ICs with Blinks:\n');

        blink=intersect ( intersect( find(SAD>=soglia_SAD),find(medie_left.*medie_right>0) ) ,...
            intersect ( find(meanK>=soglia_K),find(diff_var>0) ));

        bl_bool=1; %true if there are artifacted ICs

        if isempty(blink) %no blink component
            fprintf(file,'/ \n');
            bl_bool=0;
        else
            fprintf(file,[num2str(blink) '\n']);
            fprintf(file,'\n');    
        end

    %% Generic Discontinuities (GD)

        disp('              - Evaluating Discontinuities...')
        % ----------------------------------------------------
        % |  Writes GD header in the report file             |
        % ----------------------------------------------------

        fprintf(file,'>>>> GD - Discontinuities\n');
        fprintf(file,'Classification based on features:\n\n');
        fprintf(file,'GDSF - Generic Discontinuities Spatial Feature (threshold=%f)\n',soglia_GDSF);
        fprintf(file,'MEV - Maximum epoch variance (threshold=%f)\n\n',soglia_V);
        fprintf(file,'ICs with Discontinuities:\n');

        disc=intersect(find(GDSF>=soglia_GDSF),find(nuovaV>=soglia_V));

        dsc_bool=1; %true if there are discontinuities
        if isempty(disc) %no discontinuities
            fprintf(file,'/ \n');
            dsc_bool=0;
        else
            fprintf(file,[num2str(disc) '\n']);
            fprintf(file,'\n');    
        end

    %% Displaying results

        % ----------------------------------------------------
        % |  Write message to user: report file name         |
        % ----------------------------------------------------

        disp(' ')
        disp(['           - Results in <' out '>.'])

        fclose(file);

        %compute output variable
        art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs

        % these three are old outputs which are no more necessary in latest ADJUST version.
        soglia_D=0;
        soglia_DV=0;
        maxdin=zeros(1,size(EEG.icawinv,2));

        return

    %% The following sections have been moved to interface_ADJ in order to manage
    %% continuous data


    %% Saving artifacted ICs for further analysis

        nome=['List_' EEG.setname '.mat'];
        save (nome, 'blink', 'horiz', 'vert', 'disc');
        disp(' ')
        disp(['Artifact ICs list saved in ' nome]);


    %% IC show & remove
    % show all ICs; detected ICs are highlighted in red color. Based on
    % pop_selectcomps.

        art = nonzeros( union (union(blink,horiz) , union(vert,disc)) )'; %artifact ICs

            [EEG] = pop_selectcomps_ADJ( EEG, 1:size(EEG.icawinv,1), art, horiz, vert, blink, disc,...
                soglia_DV, diff_var, soglia_K, meanK, soglia_SED, SED, soglia_SAD, SAD, ...
                soglia_TDR, topog_DR, soglia_V, maxvar, soglia_D, maxdin );

