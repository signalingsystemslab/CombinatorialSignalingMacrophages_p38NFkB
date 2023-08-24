% Script to analyse Spearman correlations between p38 and NFkB features

% includes non-responders

%Startup 
OneDrivePath = getenv('OneDrive');
%% Load full dataset here
% AllDataSets can be found at https://doi.org/10.5281/zenodo.8274567

% load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AnalysisScripts\AllDataSets.mat'])

%% Pool the replicates
Stimuli = {'LPS', 'TNF','P3C4', 'CpG'};

for doses = 1:6
    for i = 1:2:8
        for f = 1:numel(fieldnames(DataSets.Data{i}(doses).metrics))
            FieldNames =fieldnames(DataSets.Data{i}(doses).metrics);
            DataSets.Pooled{(i+1)/2}(doses).metrics.(FieldNames{f}) = [DataSets.Data{i}(doses).metrics.(FieldNames{f});DataSets.Data{i+1}(doses).metrics.(FieldNames{f})];
        end 
    end
end
%%
% FeatureList can be found in Github repository MACKtrack - Metrics_SL here: https://github.com/signalingsystemslab/MACKtrack-for-NFkappaB-and-p38-dynamics/tree/NFkB_p38_combinatorial_signaling/Metrics_SL
NFkBFeatures= readtable([OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList'], 'Sheet', 'CorrelationsNFkB');
KTRFeatures= readtable([OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\MACKtrack_SL\Metrics_SL\FeatureList'], 'Sheet', 'CorrelationsKTR');

%% High Dose, 4 Ligands correlations between corresponding NfkB and KTR features
%Calculate correlations between all desired p38 and NFkB features 
ligand = {'LPS', 'TNF', 'P3C4', 'CpG'};
dose = 6;
n = numel(ligand); %number of conditions, eg ligands
R_value_matrix = nan(n,height(KTRFeatures));
pvalue_matrix = nan(n,height(KTRFeatures));
for m=1:1:height(KTRFeatures)
    met_nfkb= string(NFkBFeatures{m,1});
    met_ktr= string(KTRFeatures{m,1}); 
    nfkb_index= NFkBFeatures{m,2}; 
    ktr_index= KTRFeatures{m,2}; 
    count = 1;
    for i= [3,4,1,2]
        subplot(1,4, count)
        [R, pvalue] = corrplot([DataSets.Pooled{i}(dose).metrics.(met_nfkb)(:,nfkb_index), DataSets.Pooled{i}(dose).metrics.(met_ktr)(:,ktr_index)], 'varNames', {string(NFkBFeatures{m,3}), string(KTRFeatures{m,3})}, 'Type', 'Spearman','testR', 'on');
        title(ligand{i},"Interpreter", "none");
        R_value_matrix(count,m) = R(2);
        pvalue_matrix(count,m) = pvalue(2);
    count = count + 1;
    end
end


%%
figure('Position', [440,422,938,232])

% sort by mean correlation coefficient
[~, sorted_index] = sort(mean(R_value_matrix), 'descend');
heatmap(string(NFkBFeatures{sorted_index,3}),string(ligand([3,4,1,2]))',R_value_matrix(:,sorted_index))
title('Spearman correlation of NFkB and p38 features for highest dose of Pam3CSK4, CpG, LPS, TNF');
colormap([cmap([ 0 0 1],20,30,0);flipud(cmap([1 0 0],20,30,0))]); % flipud gives you two colors
caxis([-1 1]);

ax = axes;
ylim([0,4])
xlim([0,height(KTRFeatures)])
txt = '*';
for m=1:1:height(KTRFeatures)
    for i= 1:numel(ligand)
       p_value_index = sorted_index(m);
       if pvalue_matrix(i,p_value_index) <= 0.05
         text(m-0.5, i-0.4, txt, 'HorizontalAlignment','center', 'FontSize',20)
       end
    end
end
ax.Color = 'none';

set(ax, 'Ydir', 'reverse')
set(ax, 'YAxisLocation', 'Right')
% Maybe it needs some adjustment to be positioned correctly, like:
%[left bottom width height]
ax.Position(1) = 0.131;
ax.Position(2) = 0.11;
ax.Position(3) = 0.74;
ax.Position(4) = 0.8;
ax.XTick = [];
ax.YTick = [];

% saveas(gcf,[OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Correlations\HM_RVal_Spearman_Dose6_AllCells.fig'])
% saveas(gcf,[OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Correlations\HM_RVal_Spearman_Dose6_AllCells.svg'])


%% Scrambled cells control: High Dose, 4 Ligands correlations between corresponding NfkB and KTR features

ligand = {'LPS', 'TNF', 'P3C4', 'CpG'};
dose = 6;
n = numel(ligand); %number of conditions, eg ligands
R_value_matrix = nan(n,height(KTRFeatures));
pvalue_matrix = nan(n,height(KTRFeatures));
for m=1:1:height(KTRFeatures)
    figure('Position', [32,288,1404,513])
    met_nfkb= string(NFkBFeatures{m,1});
    met_ktr= string(KTRFeatures{m,1}); 
    nfkb_index= NFkBFeatures{m,2}; 
    ktr_index= KTRFeatures{m,2}; 
    count = 1;
    for i= [3,4,1,2]
        subplot(1,4, count)
        ScrambledNFkBValues = DataSets.Pooled{i}(dose).metrics.(met_nfkb)(:,nfkb_index);
        ScrambledNFkBValues = ScrambledNFkBValues(randi(length(ScrambledNFkBValues), length(ScrambledNFkBValues) ,1));
        [R, pvalue] = corrplot([ScrambledNFkBValues, DataSets.Pooled{i}(dose).metrics.(met_ktr)(:,ktr_index)], 'varNames', {string(NFkBFeatures{m,3}), string(KTRFeatures{m,3})}, 'Type', 'Spearman','testR', 'on');
        title(ligand{i},"Interpreter", "none");
        R_value_matrix_scram(count,m) = R(2);
        pvalue_matrix_scram(count,m) = pvalue(2);
    count = count + 1;
    end
    close
end

%%
figure('Position', [440,422,901,216])

% sort in order derived from anlysis with non-shuffled cell
heatmap(string(NFkBFeatures{sorted_index,3}),string(ligand([3,4,1,2]))',R_value_matrix_scram(:,sorted_index))

title('Spearman correlation of NFkB and p38 features: Neg Ctrl: Scrambled cells');
colormap([cmap([ 0 0 1],20,30,0);flipud(cmap([1 0 0],20,30,0))]); % flipud gives you two colors
caxis([-1 1]);
ax = axes;
ylim([0,4])
xlim([0,height(KTRFeatures)])
txt = '*';
for m=1:1:height(KTRFeatures)
    for i= 1:numel(ligand)
    p_value_index = sorted_index(m);
       if pvalue_matrix_scram(i,p_value_index) <= 0.05
            text(m-0.5, i-0.4, txt, 'HorizontalAlignment','center', 'FontSize',20)
       end
    end
end
ax.Color = 'none';
set(ax, 'Ydir', 'reverse')
set(ax, 'YAxisLocation', 'Right')
% Maybe it needs some adjustment to be positioned correctly, like:
%[left bottom width height]
ax.Position(1) = 0.129;
ax.Position(2) = 0.115;
ax.Position(3) = 0.74;
ax.Position(4) = 0.8;
ax.XTick = [];
ax.YTick = [];

% saveas(gcf,[OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Correlations\HM_RVal_Spearman_Dose6_AllCells_NegCtrl.fig'])
% saveas(gcf,[OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Correlations\HM_RVal_Spearman_Dose6_AllCells_NegCtrl.svg'])

%% Using 0.5 h integrals: High Dose, 4 Ligands correlations between corresponding NfkB and KTR features

window_integrals_nfkb = repmat({'intwin_p5_nfkb'},8,1);
window_integrals_ktr = repmat({'intwin_p5_ktr'},8,1);
column_index = [1:8]';
names = {'1', '2', '3', '4', '5', '6', '7','8'  };

NFkBFeaturesIntegrals = table(window_integrals_nfkb, column_index, names') ;
KTRFeaturesIntegrals = table(window_integrals_ktr, column_index, names') ;

ligand = {'LPS', 'TNF', 'P3C4', 'CpG'};
dose = 6;
n = numel(ligand); %number of conditions, eg ligands
R_value_matrix = nan(n,height(KTRFeaturesIntegrals));
pvalue_matrix = nan(n,height(KTRFeaturesIntegrals));
for m=1:1:height(KTRFeaturesIntegrals)
   figure('Position', [32,288,1404,513])
    met_nfkb= string(NFkBFeaturesIntegrals{m,1});
    met_ktr= string(KTRFeaturesIntegrals{m,1}); 
    nfkb_index= NFkBFeaturesIntegrals{m,2}; 
    ktr_index= KTRFeaturesIntegrals{m,2}; 
    count = 1;
    for i= [3,4,1,2]
        subplot(1,4, count)
        [R, pvalue] = corrplot([DataSets.Pooled{i}(dose).metrics.(met_nfkb)(:,nfkb_index), DataSets.Pooled{i}(dose).metrics.(met_ktr)(:,ktr_index)], 'varNames', {string(NFkBFeaturesIntegrals{m,3}), string(KTRFeaturesIntegrals{m,3})}, 'Type', 'Spearman','testR', 'on');
        title(ligand{i},"Interpreter", "none");
        R_value_matrix(count,m) = R(2);
        pvalue_matrix(count,m) = pvalue(2);
    count = count + 1;
    end
end

%%
figure('Position', [440,422,451,198])
heatmap(string(NFkBFeaturesIntegrals{:,3}),string(ligand([3,4,1,2]))',R_value_matrix)
title('Spearman correlation of NFkB and p38 features for highest dose of Pam3CSK4, CpG, LPS, TNF');
colormap([cmap([ 0 0 1],20,30,0);flipud(cmap([1 0 0],20,30,0))]); % flipud gives you two colors
caxis([-1 1]);

ax = axes;
ylim([0,4])
xlim([0,height(KTRFeaturesIntegrals)])
txt = '*';
for m=1:1:height(KTRFeaturesIntegrals)
    for i= 1:numel(ligand)
    p_value_index = m;
       if pvalue_matrix(i,p_value_index) <= 0.05
        text(m-0.5, i-0.4, txt, 'HorizontalAlignment','center', 'FontSize',20)
       end
    end
end
ax.Color = 'none';

set(ax, 'Ydir', 'reverse')
set(ax, 'YAxisLocation', 'Right')
% Maybe it needs some adjustment, like:
%[left bottom width height]
ax.Position(1) = 0.134;
ax.Position(2) = 0.11;
ax.Position(3) = 0.7;
ax.Position(4) = 0.785;
ax.XTick = [];
ax.YTick = [];

% saveas(gcf,[OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Correlations\HM_RVal_Spearman_Integrals_Dose6_AllCells.fig'])
% saveas(gcf,[OneDrivePath,'\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Correlations\HM_RVal_Spearman_Integrals_Dose6_AllCells.svg'])
