%% Run Script + function to calculate mutual information from p38 only, NFkB only, or p38 and p38+NFKB dynamic features

% "Mutual information (MI) was calculated from p38, NFκB, or p38 and NFκB dynamic features 
% from 2 pooled biological replicates from the indicated sets of stimulation conditions. 
% Prior to calculation, dynamic feature containing NaN values were removed and 
% the features were z-scored. MI was extrapolated using jackknife resampling to control for different sample sizes 
% with 12 sample subsets, each containing 65%-90% of total samples (Adelaja et al, 2021). 
% MI was calculated using a published method optimized for determining mutual information 
% between discrete and continuous data sets (Ross, 2014)."

% Using input data tables found HERE
% in this repository:
% 'CombinatorialSignalingMacrophages_p38NFkB\MutualInformation\Inputs'
% 'doses' indicate position in input data structure for respective stimulation condition (ligand+dose)
%% 4 Stimuli, Dose 6
doses = [6, 12, 18, 24];
[AllStim_PooledRep_Dose6_Ross(1).info_extrap,AllStim_PooledRep_Dose6_Ross(1).info_fullset,AllStim_PooledRep_Dose6_Ross(1).info_subsets,AllStim_PooledRep_Dose6_Ross(1).weights,  AllStim_PooledRep_Dose6_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[AllStim_PooledRep_Dose6_Ross(2).info_extrap,AllStim_PooledRep_Dose6_Ross(2).info_fullset,AllStim_PooledRep_Dose6_Ross(2).info_subsets,AllStim_PooledRep_Dose6_Ross(2).weights,  AllStim_PooledRep_Dose6_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[AllStim_PooledRep_Dose6_Ross(3).info_extrap,AllStim_PooledRep_Dose6_Ross(3).info_fullset,AllStim_PooledRep_Dose6_Ross(3).info_subsets,AllStim_PooledRep_Dose6_Ross(3).weights,  AllStim_PooledRep_Dose6_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

%% 3 PAMPs, Dose 6
doses = [6, 18, 24];
[PAMPs_PooledRep_Dose6_Ross(1).info_extrap,PAMPs_PooledRep_Dose6_Ross(1).info_fullset,PAMPs_PooledRep_Dose6_Ross(1).info_subsets,PAMPs_PooledRep_Dose6_Ross(1).weights,  PAMPs_PooledRep_Dose6_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[PAMPs_PooledRep_Dose6_Ross(2).info_extrap,PAMPs_PooledRep_Dose6_Ross(2).info_fullset,PAMPs_PooledRep_Dose6_Ross(2).info_subsets,PAMPs_PooledRep_Dose6_Ross(2).weights,  PAMPs_PooledRep_Dose6_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[PAMPs_PooledRep_Dose6_Ross(3).info_extrap,PAMPs_PooledRep_Dose6_Ross(3).info_fullset,PAMPs_PooledRep_Dose6_Ross(3).info_subsets,PAMPs_PooledRep_Dose6_Ross(3).weights,  PAMPs_PooledRep_Dose6_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

%% LPS Mock vs Dose 6
doses = [1,6];
[LPS_PooledRep_MOCKvsDose6_Ross(1).info_extrap,LPS_PooledRep_MOCKvsDose6_Ross(1).info_fullset,LPS_PooledRep_MOCKvsDose6_Ross(1).info_subsets,LPS_PooledRep_MOCKvsDose6_Ross(1).weights,  LPS_PooledRep_MOCKvsDose6_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[LPS_PooledRep_MOCKvsDose6_Ross(2).info_extrap,LPS_PooledRep_MOCKvsDose6_Ross(2).info_fullset,LPS_PooledRep_MOCKvsDose6_Ross(2).info_subsets,LPS_PooledRep_MOCKvsDose6_Ross(2).weights,  LPS_PooledRep_MOCKvsDose6_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[LPS_PooledRep_MOCKvsDose6_Ross(3).info_extrap,LPS_PooledRep_MOCKvsDose6_Ross(3).info_fullset,LPS_PooledRep_MOCKvsDose6_Ross(3).info_subsets,LPS_PooledRep_MOCKvsDose6_Ross(3).weights,  LPS_PooledRep_MOCKvsDose6_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);
%% TNF Mock vs Dose 6
doses = [7,12];
[TNF_PooledRep_MOCKvsDose6_Ross(1).info_extrap,TNF_PooledRep_MOCKvsDose6_Ross(1).info_fullset,TNF_PooledRep_MOCKvsDose6_Ross(1).info_subsets,TNF_PooledRep_MOCKvsDose6_Ross(1).weights,  TNF_PooledRep_MOCKvsDose6_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[TNF_PooledRep_MOCKvsDose6_Ross(2).info_extrap,TNF_PooledRep_MOCKvsDose6_Ross(2).info_fullset,TNF_PooledRep_MOCKvsDose6_Ross(2).info_subsets,TNF_PooledRep_MOCKvsDose6_Ross(2).weights,  TNF_PooledRep_MOCKvsDose6_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[TNF_PooledRep_MOCKvsDose6_Ross(3).info_extrap,TNF_PooledRep_MOCKvsDose6_Ross(3).info_fullset,TNF_PooledRep_MOCKvsDose6_Ross(3).info_subsets,TNF_PooledRep_MOCKvsDose6_Ross(3).weights,  TNF_PooledRep_MOCKvsDose6_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);
%% P3C4 Mock vs Dose 6
% 20230322: rerun iwth 12 subsets, saved
doses = [13,18];
[P3C4_PooledRep_MOCKvsDose6_Ross(1).info_extrap,P3C4_PooledRep_MOCKvsDose6_Ross(1).info_fullset,P3C4_PooledRep_MOCKvsDose6_Ross(1).info_subsets,P3C4_PooledRep_MOCKvsDose6_Ross(1).weights,  P3C4_PooledRep_MOCKvsDose6_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[P3C4_PooledRep_MOCKvsDose6_Ross(2).info_extrap,P3C4_PooledRep_MOCKvsDose6_Ross(2).info_fullset,P3C4_PooledRep_MOCKvsDose6_Ross(2).info_subsets,P3C4_PooledRep_MOCKvsDose6_Ross(2).weights,  P3C4_PooledRep_MOCKvsDose6_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[P3C4_PooledRep_MOCKvsDose6_Ross(3).info_extrap,P3C4_PooledRep_MOCKvsDose6_Ross(3).info_fullset,P3C4_PooledRep_MOCKvsDose6_Ross(3).info_subsets,P3C4_PooledRep_MOCKvsDose6_Ross(3).weights,  P3C4_PooledRep_MOCKvsDose6_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);
%% CpG Mock vs Dose 6
doses = [19,24];
[CpG_PooledRep_MOCKvsDose6_Ross(1).info_extrap,CpG_PooledRep_MOCKvsDose6_Ross(1).info_fullset,CpG_PooledRep_MOCKvsDose6_Ross(1).info_subsets,CpG_PooledRep_MOCKvsDose6_Ross(1).weights,  CpG_PooledRep_MOCKvsDose6_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[CpG_PooledRep_MOCKvsDose6_Ross(2).info_extrap,CpG_PooledRep_MOCKvsDose6_Ross(2).info_fullset,CpG_PooledRep_MOCKvsDose6_Ross(2).info_subsets,CpG_PooledRep_MOCKvsDose6_Ross(2).weights,  CpG_PooledRep_MOCKvsDose6_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[CpG_PooledRep_MOCKvsDose6_Ross(3).info_extrap,CpG_PooledRep_MOCKvsDose6_Ross(3).info_fullset,CpG_PooledRep_MOCKvsDose6_Ross(3).info_subsets,CpG_PooledRep_MOCKvsDose6_Ross(3).weights,  CpG_PooledRep_MOCKvsDose6_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

%% 4x Ligand dose responses incl mock, 

%LPS Doses 1:6
doses = 1:6;
[LPS_PooledRep_AllDoses_inclMock_Ross(1).info_extrap,LPS_PooledRep_AllDoses_inclMock_Ross(1).info_fullset,LPS_PooledRep_AllDoses_inclMock_Ross(1).info_subsets,LPS_PooledRep_AllDoses_inclMock_Ross(1).weights,  LPS_PooledRep_AllDoses_inclMock_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[LPS_PooledRep_AllDoses_inclMock_Ross(2).info_extrap,LPS_PooledRep_AllDoses_inclMock_Ross(2).info_fullset,LPS_PooledRep_AllDoses_inclMock_Ross(2).info_subsets,LPS_PooledRep_AllDoses_inclMock_Ross(2).weights,  LPS_PooledRep_AllDoses_inclMock_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[LPS_PooledRep_AllDoses_inclMock_Ross(3).info_extrap,LPS_PooledRep_AllDoses_inclMock_Ross(3).info_fullset,LPS_PooledRep_AllDoses_inclMock_Ross(3).info_subsets,LPS_PooledRep_AllDoses_inclMock_Ross(3).weights,  LPS_PooledRep_AllDoses_inclMock_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

% TNF Doses 1:6 
doses = 7:12;
[TNF_PooledRep_AllDoses_inclMock_Ross(1).info_extrap,TNF_PooledRep_AllDoses_inclMock_Ross(1).info_fullset,TNF_PooledRep_AllDoses_inclMock_Ross(1).info_subsets,TNF_PooledRep_AllDoses_inclMock_Ross(1).weights,  TNF_PooledRep_AllDoses_inclMock_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[TNF_PooledRep_AllDoses_inclMock_Ross(2).info_extrap,TNF_PooledRep_AllDoses_inclMock_Ross(2).info_fullset,TNF_PooledRep_AllDoses_inclMock_Ross(2).info_subsets,TNF_PooledRep_AllDoses_inclMock_Ross(2).weights,  TNF_PooledRep_AllDoses_inclMock_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[TNF_PooledRep_AllDoses_inclMock_Ross(3).info_extrap,TNF_PooledRep_AllDoses_inclMock_Ross(3).info_fullset,TNF_PooledRep_AllDoses_inclMock_Ross(3).info_subsets,TNF_PooledRep_AllDoses_inclMock_Ross(3).weights,  TNF_PooledRep_AllDoses_inclMock_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

% P3C4 Doses 1:6
doses = 13:18;
[P3C4_PooledRep_AllDoses_inclMock_Ross(1).info_extrap,P3C4_PooledRep_AllDoses_inclMock_Ross(1).info_fullset,P3C4_PooledRep_AllDoses_inclMock_Ross(1).info_subsets,P3C4_PooledRep_AllDoses_inclMock_Ross(1).weights,  P3C4_PooledRep_AllDoses_inclMock_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[P3C4_PooledRep_AllDoses_inclMock_Ross(2).info_extrap,P3C4_PooledRep_AllDoses_inclMock_Ross(2).info_fullset,P3C4_PooledRep_AllDoses_inclMock_Ross(2).info_subsets,P3C4_PooledRep_AllDoses_inclMock_Ross(2).weights,  P3C4_PooledRep_AllDoses_inclMock_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[P3C4_PooledRep_AllDoses_inclMock_Ross(3).info_extrap,P3C4_PooledRep_AllDoses_inclMock_Ross(3).info_fullset,P3C4_PooledRep_AllDoses_inclMock_Ross(3).info_subsets,P3C4_PooledRep_AllDoses_inclMock_Ross(3).weights,  P3C4_PooledRep_AllDoses_inclMock_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

% CpG Doses 1:6
doses = 19:24;
[CpG_PooledRep_AllDoses_inclMock_Ross(1).info_extrap,CpG_PooledRep_AllDoses_inclMock_Ross(1).info_fullset,CpG_PooledRep_AllDoses_inclMock_Ross(1).info_subsets,CpG_PooledRep_AllDoses_inclMock_Ross(1).weights,  CpG_PooledRep_AllDoses_inclMock_Ross(1).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_both', 'CC_input_AllStim_Rep2_AllDoses_both', 'nfkb', doses);
[CpG_PooledRep_AllDoses_inclMock_Ross(2).info_extrap,CpG_PooledRep_AllDoses_inclMock_Ross(2).info_fullset,CpG_PooledRep_AllDoses_inclMock_Ross(2).info_subsets,CpG_PooledRep_AllDoses_inclMock_Ross(2).weights,  CpG_PooledRep_AllDoses_inclMock_Ross(2).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_NFkB', 'CC_input_AllStim_Rep2_AllDoses_NFkB', 'nfkb', doses);
[CpG_PooledRep_AllDoses_inclMock_Ross(3).info_extrap,CpG_PooledRep_AllDoses_inclMock_Ross(3).info_fullset,CpG_PooledRep_AllDoses_inclMock_Ross(3).info_subsets,CpG_PooledRep_AllDoses_inclMock_Ross(3).weights,  CpG_PooledRep_AllDoses_inclMock_Ross(3).DataFrame] = get_mutual_info_Ross('CC_input_AllStim_Rep1_AllDoses_KTR','CC_input_AllStim_Rep2_AllDoses_KTR', 'ktr', doses);

%% Function that calls jacknife function (which calls mututal information calculation function) for respective inputs
function [fitI,info1,I,Q, DataFrame] = get_mutual_info_Ross(DataSetName_Rep1,DataSetName_Rep2, reporter, doses)
OneDrivePath = getenv('OneDrive');

% Read dynamic features for Replicate 1
% adjust load path as needed
all_metrics1 = load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Information theory\inforun_inputs\', DataSetName_Rep1, '.mat']);
all_metrics1 = all_metrics1.nfkb;

DataFrame1 = table();
for i = 1:length(all_metrics1)
    DataFramePart = struct2table(all_metrics1(i).metrics);
    if strcmp(reporter ,'nfkb')
        DataFramePart(:,end+1) = repmat({i}, length(all_metrics1(i).metrics.responder_status_nfkb1),1);
    else       
        DataFramePart(:,end+1) = repmat({i}, length(all_metrics1(i).metrics.responder_status_ktr1),1);
    end
    DataFramePart = [DataFramePart(:,end),DataFramePart(:,1:end-1)];
    DataFrame1 = [DataFrame1; DataFramePart];
end

% Read dynamic features for Replicate 2
% adjust load path as needed
all_metrics2 = load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Information theory\inforun_inputs\', DataSetName_Rep2, '.mat']);
all_metrics2 = all_metrics2.nfkb;

DataFrame2 = table();
for i = 1:length(all_metrics2)
    DataFramePart = struct2table(all_metrics2(i).metrics);
    if strcmp(reporter ,'nfkb')
        DataFramePart(:,end+1) = repmat({i}, length(all_metrics2(i).metrics.responder_status_nfkb1),1);
    else       
        DataFramePart(:,end+1) = repmat({i}, length(all_metrics2(i).metrics.responder_status_ktr1),1);
    end
    DataFramePart = [DataFramePart(:,end),DataFramePart(:,1:end-1)];
    DataFrame2 = [DataFrame2; DataFramePart];
end

% Pool the 2 biological replicates
DataFrame = [DataFrame1;DataFrame2];

DataFrame = table2array(DataFrame);
DataFrame = DataFrame(ismember(DataFrame(:,1), doses),:);

% Remove dynamic features with NaN
DataFrame= DataFrame(:,all(~isnan(DataFrame))); 

% Z score dynamic features
DataFrame = [DataFrame(:,1), zscore(DataFrame(:,2:end))];

[fitI,info1,I,Q] = jacknife_forRoss2014MI(DataFrame(1:end,1)', DataFrame(1:end,2:end)',9,12, 'verbose');
 
% adjust save path as needed to save output structure
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Information theory\Ross_2014\Jackknife_diagnosis\PooledReps\',DataSetName_Rep1, '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Information theory\Ross_2014\Jackknife_diagnosis\PooledReps\',DataSetName_Rep1,'.svg'])
end