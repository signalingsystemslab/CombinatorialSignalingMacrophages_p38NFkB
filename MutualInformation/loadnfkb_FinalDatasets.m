function [nfkb, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,filename, varargin)

%% This function generates one large input structure from p38 / NFkB dynamic features of multiple experimental conditions for input to mutual information calculation

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% NOTE: originally written to generate data sets for running information theory calculations in Adelaja, Taylor et al., 2021, Immunity, with automatic loading of AllMeasurements.mat output from Macktrac
% [nfkb, all_dims, names_1D] = loadnfkb()
% LOADNFKB loads the master set of NFkB runs used for (most) BMDM analyses
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

% In this project:
% Using input data set ("AllDataSets.mat") found at https://doi.org/10.5281/zenodo.8274567

% Example of how to call this function for the p38 NFkB combinatorial signaling project:
    % "DataSets" is loaded from "AllDataSets.mat"
    % stimuli = {'LPS', 'TNF', 'P3C4','CpG'};
    % For each input file to be generated, this function is called using one of the commands like the examples below
    % [CC_input_AllStim_Rep1_AllDoses_KTR, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,'CC_input_AllStim_Rep1_AllDoses_KTR','KTRonly', 1, 'NFkBonly', 0, 'Stimuli', stimuli,'doses', [1:6]); 
    % [CC_input_AllStim_Rep1_AllDoses_NFkB, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,'CC_input_AllStim_Rep1_AllDoses_NFkB','KTRonly', 0, 'NFkBonly', 1, 'Stimuli', stimuli,'doses', [1:6]);
    % [CC_input_AllStim_Rep1_AllDoses_both, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,'CC_input_AllStim_Rep1_AllDoses_both','KTRonly', 0, 'NFkBonly', 0, 'Stimuli', stimuli,'doses', [1:6]);
    % [CC_input_AllStim_Rep2_AllDoses_KTR, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,'CC_input_AllStim_Rep2_AllDoses_KTR','KTRonly', 1, 'NFkBonly', 0, 'Stimuli', stimuli,'doses', [1:6]); 
    % [CC_input_AllStim_Rep2_AllDoses_NFkB, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,'CC_input_AllStim_Rep2_AllDoses_NFkB','KTRonly', 0, 'NFkBonly', 1, 'Stimuli', stimuli,'doses', [1:6]);
    % [CC_input_AllStim_Rep2_AllDoses_both, all_dims, names_1D] = loadnfkb_FinalDatasets(DataSets,'CC_input_AllStim_Rep2_AllDoses_both','KTRonly', 0, 'NFkBonly', 0, 'Stimuli', stimuli,'doses', [1:6]); 
% -----------------------------------------------------------------------------------------------------------------------

p=inputParser; 
addParameter(p, 'NFkBonly', 0)
addParameter(p, 'KTRonly', 0)
addParameter(p, 'Stimuli', {'LPS', 'TNF', 'P3C4', 'CpG'})
addParameter(p, 'doses', [1:6])

parse(p,varargin{:})

NFkBonly = p.Results.NFkBonly;
KTRonly = p.Results.KTRonly;
%%
% Define name of nfkb file
P = mfilename('fullpath');
P2 = mfilename;
nfkb_name = [P(1:(length(P)-length(P2))) , filename,'.mat'];

%%  
%Startup 
OneDrivePath = getenv('OneDrive');

stimuli = p.Results.Stimuli;
doses = p.Results.doses;

% Combine processed data into a single structure
if ~exist(nfkb_name,'file')

% NOTE: Edit load path as needed to load input data structure
% Using input data set ("AllDataSets.mat") found at https://doi.org/10.5281/zenodo.8274567
%    load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AnalysisScripts\AllDataSets.mat'])

    % Load list of features to be included in calculations from appropriate path
    % FeatureList.xlsx can be found in MACKtrack here: https://github.com/signalingsystemslab/MACKtrack-for-NFkappaB-and-p38-dynamics/tree/NFkB_p38_combinatorial_signaling/Metrics_SL
    feature_table_nfkb= readtable([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList.xlsx'], 'Sheet', 'AllFeatNFkB_cc');
    feature_table_ktr= readtable([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList.xlsx'], 'Sheet', 'AllFeatKTR_cc');

    features_nfkb= table2cell(feature_table_nfkb(:,1)); %list of metrics/features to be plotted
    feature_index_nfkb = table2cell(feature_table_nfkb(:,2)); %index for metrics with multiple columns (eg duration, etc.)
    features_ktr= table2cell(feature_table_ktr(:,1)); %list of metrics/features to be plotted
    feature_index_ktr= table2cell(feature_table_ktr(:,2)); %index for metrics with multiple columns (eg duration, etc.)
    nfkb = struct;
    ids = cell(1, 1);
    for f = 1:numel(features_nfkb)
        count = 1;
        for s = 1:numel(stimuli)
            switch stimuli{s}
 % NOTE: TOGGLE this to generate structure for replicate 1 or 2
                case 'LPS'
%                    RepNames = {'LPS_R1', 'LPS_R2'};
%                    RepNames = {'LPS_R1'};
                    RepNames = {'LPS_R2'};
                case 'TNF'
%                    RepNames = {'TNF_R1', 'TNF_R2'};
%                    RepNames = {'TNF_R1'};
                    RepNames = {'TNF_R2'};
                case 'P3C4'
%                    RepNames = {'P3C4_R1', 'P3C4_R2'};
 %                    RepNames = {'P3C4_R1'};
                     RepNames = {'P3C4_R2'};
                case 'CpG'
%                  RepNames = {'CpG_R1'};
                  RepNames = {'CpG_R2'};
%                    RepNames = {'CpG_R1', 'CpG_R2'};
            end
            for r = 1:numel(RepNames)
                RepIndex = find(contains(DataSets.RepNames, RepNames{r}));
                for d = 1:numel(doses)
                    nfkb(count).metrics.([features_nfkb{f},num2str(feature_index_nfkb{f})]) = DataSets.Data{RepIndex}(doses(d)).metrics.(features_nfkb{f})(:,feature_index_nfkb{f});            
                    nfkb(count).metrics.([features_ktr{f},num2str(feature_index_ktr{f})]) = DataSets.Data{RepIndex}(doses(d)).metrics.(features_ktr{f})(:,feature_index_ktr{f});  
                    ids{count} = [RepNames{r},'_', num2str(doses(d))]; %strcat(stimuli{i},doses{j});
                    nfkb(count).id = ids(count);
                    count = count + 1;

                end
            end
        end
    end
    %include this to remove nfkb/ktr to test only ktr/nfkb
    if NFkBonly == 1
        names = fieldnames(nfkb(1).metrics);
        for i =1:length(ids)
            nfkb(i).metrics = rmfield(nfkb(i).metrics,names(contains(names, 'ktr'))) ;
        end
    end
    
    if KTRonly == 1
        names = fieldnames(nfkb(1).metrics);
        for i =1:length(ids)
            nfkb(i).metrics = rmfield(nfkb(i).metrics,names(contains(names, 'nfkb'))) ;
        end
    end

    save(nfkb_name, 'nfkb')
else
    load(nfkb_name);
end


% If extra output arguments are defined, combine all metrics/measurements into one big matrix 
% (record dimension names as well)
if nargout~=1
    mfields = fieldnames(nfkb(1).metrics);
    names_1D = {};
    for i = 1:length(nfkb)
        % Combine all dimensions into a single gigantic cell thing
        all_dims{i} = [];
        for j = 1:length(mfields)
            all_dims{i} = cat(2,all_dims{i},nfkb(i).metrics.(mfields{j}));
            if i==1
                for k = 1:size(nfkb(i).metrics.(mfields{j}),2)
%                    names_1D = cat(1,names_1D,[mfields{j},'_',numseq(k,3)]);
                    names_1D = cat(1,names_1D,[mfields{j}]);
                end
            end
        end
    end

    % Reload full nfkb.mat
    load(nfkb_name)
end

% If called without output, assign everything into base workspace
if nargout<1
    assignin('base','nfkb',nfkb)
    assignin('base','all_dims',all_dims)
    assignin('base','names_1D',names_1D)
end