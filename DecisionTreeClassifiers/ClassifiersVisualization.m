% Script to visualize 5-fold cross validation results of Decision Tree
% Ensemble classification of p38 and NFkB dynamic features by ligand
% identity and dose

%Startup 
OneDrivePath = getenv('OneDrive');

%% Loading Models with same numbers of cells per class
% Trained models with evaluation scores can be found here:  doi.org/10.5281/zenodo.8274744
% load as appropriate
%load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z score\PooledReps\ModelResults_Both_PooledReps_EqualCell.mat']);
%load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z score\PooledReps\ModelResults_NFkB_PooledReps_EqualCell.mat']);
%load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z score\PooledReps\ModelResults_KTR_PooledReps_EqualCell.mat']);
%load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z score\PooledReps\ModelResults_BothNFkBscram_PooledReps_EqualCell.mat']);
%load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z score\PooledReps\ModelResults_BothKTRscram_PooledReps_EqualCell.mat']);

%% Figure 2 AllStimDose6 p38only Confusion Matrix

name =  'AllStimHighDosePool';
order = ["3", "4","1","2"]; %P3C4, CpG, LPS, TNF
figure('Position',[391, 156,481, 420])
Val_Matrix = confusionchart(ModelKTR.(name).ClassificationEnsemble.Y,validationPredictionsKTR.(name)  );
Val_Matrix.RowSummary = 'row-normalized';
Val_Matrix.ColumnSummary = 'column-normalized';
Val_Matrix.Normalization= 'total-normalized';
sortClasses(Val_Matrix, order)
title([name, ', 5-fold cross-val'])
saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\KTR_',name,'_5fVal_EqualCell' '.fig'])
saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\KTR_',name,'_5fVal_EqualCell' '.svg'])
%% Figure 2 AllStimDose6 p38only Validation Accuracy
name =  'AllStimHighDosePool';
validationAccuracyKTR.(name)
%% Figure 2 AllStimDose6 p38only F1 Scores by class
name =  'AllStimHighDosePool';
score_names = 'Fscore';
order = [3, 4,1,2]; %P3C4, CpG, LPS, TNF

FieldNames =fieldnames(ModelKTR);
        plot_data_ktr = StatsValidationKTR.(name).(score_names)(order)';
        plot_data = plot_data_ktr;
        
        figure('Position',[391, 156,240, 227]);
        b= bar(plot_data, 0.75, 'FaceColor','flat');
        ylabel('F1-score')
        xlabel('Stimulus')
        xticklabels(["P3C4", "CpG", "LPS", "TNF"])
        b(1).FaceColor = [0.266666666666667	0.447058823529412	0.768627450980392];
        ylim([0 1])
        xlim([0.3 4.7])
        title(['5-fold cross-validation'])

        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels1 = string(round(b(1).YData,2));
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')

%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\KTR_',name,'_F1score_Class' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\KTR_',name,'_F1score_Class' '.svg'])
%% Figure 3 AllStimDose6 p38only, NFkB only, both, scrambled controls  Confusion Matrix
name =  'AllStimHighDosePool';
order = ["3", "4","1","2"]; %P3C4, CpG, LPS, TNF
figure('Position',[391, 156,2000, 444])
t =tiledlayout(1,5);
   for j = 1:5
    nexttile
   switch j
       case 1
       Val_Matrix = confusionchart(ModelKTR.(name).ClassificationEnsemble.Y,validationPredictionsKTR.(name)  );
       case 2
       Val_Matrix = confusionchart(ModelNFkB.(name).ClassificationEnsemble.Y,validationPredictionsNFkB.(name)  );
       case 3
       Val_Matrix = confusionchart(ModelBoth.(name).ClassificationEnsemble.Y,validationPredictionsBoth.(name)  );
       case 4
       Val_Matrix = confusionchart(ModelBothNFkBscram.(name).ClassificationEnsemble.Y,validationPredictionsBothNFkBscram.(name)  );
       case 5
       Val_Matrix = confusionchart(ModelBothKTRscram.(name).ClassificationEnsemble.Y,validationPredictionsBothKTRscram.(name)  );
   end
    Val_Matrix.RowSummary = 'row-normalized';
    Val_Matrix.ColumnSummary = 'column-normalized';
    Val_Matrix.Normalization= 'total-normalized';
    sortClasses(Val_Matrix, order)
   end
title(t,['Feature-based classification', ', 5-fold cross-val'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_5fVal_EqualCell' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_5fVal_EqualCell' '.svg'])

%% Figure 3 AllStimDose6 p38only, NFkB only, both, both scrambled  Validation Accuracy, incl graph
name =  'AllStimHighDosePool';
   for j = 1:5
   switch j
       case 1
        validationAccuracy(j)= validationAccuracyKTR.(name);
       case 2
        validationAccuracy(j)= validationAccuracyNFkB.(name);
       case 3
        validationAccuracy(j)= validationAccuracyBoth.(name);
       case 4
        validationAccuracy(j)= validationAccuracyBothNFkBscram.(name);
       case 5
        validationAccuracy(j)= validationAccuracyBothKTRscram.(name);
   end
   end
   validationAccuracy

   % plot data
   colors = setcolors;
   plot_data = validationAccuracy;
        figure('Position',[391, 156,267, 264]);
        x=1;
        b= bar(x,plot_data', 0.75, 'FaceColor','flat');
        ylabel('Overall Accuracy')
        xlabel('Stimulus')
       b(1).FaceColor = colors.p38;
       b(2).FaceColor = colors.mVen;
       b(3).FaceColor = colors.mVenKTR;
       b(4).FaceColor = colors.NFkBscram;
       b(5).FaceColor = colors.p38scram;
        ylim([0.5 1])
        title(['Model accuracy'])

    for i = 1:5
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(round(b(i).YData,2));
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    end
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All5_', name, 'ValidationAccuracies', '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All5_',name, 'ValidationAccuracies','.svg'])


   %% Figure 3 AllStimDose6 p38, NFkB only, Both, both scram F1 Scores by class (incl TNF only)
name =  'AllStimHighDosePool';
score_names = 'Fscore';
order = [3, 4,1,2]; %P3C4, CpG, LPS, TNF
colors = setcolors;

 plot_data_KTR = StatsValidationKTR.(name).(score_names)(order)';
 plot_data_nfkb = StatsValidationNFkB.(name).(score_names)(order)';
 plot_data_both = StatsValidationBoth.(name).(score_names)(order)';
 plot_data_bothNFkBscram = StatsValidationBothNFkBscram.(name).(score_names)(order)';
 plot_data_bothKTRscram = StatsValidationBothKTRscram.(name).(score_names)(order)';
 plot_data = [plot_data_KTR, plot_data_nfkb, plot_data_both, plot_data_bothNFkBscram, plot_data_bothKTRscram];  
 plot_data = reshape(plot_data, [],5);

 figure('Position',[391, 156,267, 264]);
% TNF class only
        x =1;
        b= bar(x,plot_data(4,:), 0.75, 'FaceColor','flat');
        ylabel('F1-score')
        xlabel('Stimulus')
        xticklabels(["TNF"])
        b(1).FaceColor = colors.p38;
        b(2).FaceColor = colors.mVen;
        b(3).FaceColor = colors.mVenKTR;
        b(4).FaceColor = colors.NFkBscram;
        b(5).FaceColor = colors.p38scram;
         
        ylim([0 1])
        ylim([0.5 1])
        title(['5-fold cross-validation'])

 for i = 1:5
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(round(b(i).YData,2));
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
 end

%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_F1score_Class_TNF' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_F1score_Class_TNF' '.svg'])
%% Figure 4 PAMPs p38 only, NFkB only, Both, both scrambled F1 Scores by class
name =  'PAMPsHighDosePool';
score_names = 'Fscore';
order = [2, 3,1]; %P3C4, CpG, LPS
colors = setcolors;

 plot_data_ktr = StatsValidationKTR.(name).(score_names)(order)';
 plot_data_nfkb = StatsValidationNFkB.(name).(score_names)(order)';
 plot_data_both = StatsValidationBoth.(name).(score_names)(order)';
 plot_data_bothNFkBscram = StatsValidationBothNFkBscram.(name).(score_names)(order)';
 plot_data_bothKTRscram = StatsValidationBothKTRscram.(name).(score_names)(order)';
 plot_data = [plot_data_ktr, plot_data_nfkb, plot_data_both, plot_data_bothNFkBscram, plot_data_bothKTRscram];  
 plot_data = reshape(plot_data, [],5);

        figure('Position',[391, 156,760, 365]);
        b= bar(plot_data, 0.75, 'FaceColor','flat');
        ylabel('F1-score')
        xlabel('Stimulus')
        xticklabels(["P3C4", "CpG", "LPS"])
        b(1).FaceColor = colors.p38;
        b(2).FaceColor = colors.mVen;
        b(3).FaceColor = colors.mVenKTR;
        b(4).FaceColor = colors.NFkBscram;
        b(5).FaceColor = colors.p38scram;
         
        ylim([0.5 1])
        title(['5-fold cross-validation'])

 for i = 1:5
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(round(b(i).YData,2));
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
 end

%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_F1score_Class' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_F1score_Class' '.svg'])
%% Figure 4 PAMPs Dose6 p38only, NFkB only, both  Confusion Matrix

name =  'PAMPsHighDosePool';
order = ["3", "4","1"]; %P3C4, CpG, LPS
figure('Position',[391, 156,1924, 372])
t =tiledlayout(1,5);
   for j = 1:5
    nexttile
   switch j
       case 1
       Val_Matrix = confusionchart(ModelKTR.(name).ClassificationEnsemble.Y,validationPredictionsKTR.(name)  );
       case 2
       Val_Matrix = confusionchart(ModelNFkB.(name).ClassificationEnsemble.Y,validationPredictionsNFkB.(name)  );
       case 3
       Val_Matrix = confusionchart(ModelBoth.(name).ClassificationEnsemble.Y,validationPredictionsBoth.(name)  );
       case 4
       Val_Matrix = confusionchart(ModelBothNFkBscram.(name).ClassificationEnsemble.Y,validationPredictionsBothNFkBscram.(name)  );
       case 5
       Val_Matrix = confusionchart(ModelBothKTRscram.(name).ClassificationEnsemble.Y,validationPredictionsBothKTRscram.(name)  );
   end
    Val_Matrix.RowSummary = 'row-normalized';
    Val_Matrix.ColumnSummary = 'column-normalized';
    Val_Matrix.Normalization= 'total-normalized';
    sortClasses(Val_Matrix, order)
   end
title(t,['Feature-based classification', ', 5-fold cross-val'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_5fVal_EqualCell' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All5_',name,'_5fVal_EqualCell' '.svg'])

%% Figure 4 PAMPs p38only, NFkB only, both, both srambled  Validation Accuracy, incl graph
name =  'PAMPsHighDosePool';
  for j = 1:5
   switch j
       case 1
        validationAccuracy(j)= validationAccuracyKTR.(name);
       case 2
        validationAccuracy(j)= validationAccuracyNFkB.(name);
       case 3
        validationAccuracy(j)= validationAccuracyBoth.(name);
       case 4
        validationAccuracy(j)= validationAccuracyBothNFkBscram.(name);
       case 5
        validationAccuracy(j)= validationAccuracyBothKTRscram.(name);
   end
   end
   validationAccuracy

   % plot data
   colors = setcolors;
   plot_data = validationAccuracy;
        figure('Position',[391, 156,267, 264]);
        x=1;
        b= bar(x,plot_data', 0.75, 'FaceColor','flat');
        ylabel('Overall Accuracy')
        xlabel('Stimulus')
       b(1).FaceColor = colors.p38;
       b(2).FaceColor = colors.mVen;
       b(3).FaceColor = colors.mVenKTR;
       b(4).FaceColor = colors.NFkBscram;
       b(5).FaceColor = colors.p38scram;
       ylim([0.5 1])
       title(['Model accuracy'])

    for i = 1:5
        xtips = b(i).XEndPoints;
        ytips = b(i).YEndPoints;
        labels = string(round(b(i).YData,2));
        text(xtips,ytips,labels,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    end
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All5_', name, 'ValidationAccuracies', '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All5_',name, 'ValidationAccuracies','.svg'])
%% Figure 5 Dose Responses Dose1-6 p38only, NFkB only, both  Confusion Matrix
name =  {'P3C4_Doses_Pool','CpG_Doses_Pool','LPS_Doses_Pool','TNF_Doses_Pool'};
for i = 1:numel(name)
    figure('Position',[391, 156,1520, 444])
    t =tiledlayout(1,3);
       for j = 1:3
        nexttile
       switch j
           case 1
           Val_Matrix = confusionchart(ModelKTR.(name{i}).ClassificationEnsemble.Y,validationPredictionsKTR.(name{i})  );
           case 2
           Val_Matrix = confusionchart(ModelNFkB.(name{i}).ClassificationEnsemble.Y,validationPredictionsNFkB.(name{i})  );
           case 3
           Val_Matrix = confusionchart(ModelBoth.(name{i}).ClassificationEnsemble.Y,validationPredictionsBoth.(name{i})  );
       end
        Val_Matrix.RowSummary = 'row-normalized';
        Val_Matrix.ColumnSummary = 'column-normalized';
        Val_Matrix.Normalization= 'total-normalized';
       end
    title(t,['Feature-based classification', ', 5-fold cross-val ', name{i}])
%    saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\All3_',name{i},'_5fVal_EqualCell' '.fig'])
%    saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\All3_',name{i},'_5fVal_EqualCell' '.svg'])
end

%% Figure 5 Dose responses p38only, NFkB only, both  Validation Accuracy, including bar graph
name =  {'P3C4_Doses_Pool','CpG_Doses_Pool','LPS_Doses_Pool','TNF_Doses_Pool'};
for i = 1:numel(name)
   for j = 1:3
   switch j
       case 1
        validationAccuracy(i,j)= validationAccuracyKTR.(name{i});
       case 2
        validationAccuracy(i,j)= validationAccuracyNFkB.(name{i});
        case 3
        validationAccuracy(i,j)= validationAccuracyBoth.(name{i});
   end
   end
end
   validationAccuracy
% plot data
   colors = setcolors;
   plot_data = validationAccuracy;
        figure('Position',[391, 156,762, 227]);
        b= bar(plot_data, 0.75, 'FaceColor','flat');
        ylabel('Overall Accuracy')
        xlabel('Stimulus')
        xticklabels(["P3C4", "CpG", "LPS", "TNF"])
        b(1).FaceColor = colors.p38;
        b(2).FaceColor = colors.mVen;
        b(3).FaceColor = colors.mVenKTR;  

        ylim([0 1])
        title(['5-fold cross-validation'])

        xtips1 = b(1).XEndPoints;
        ytips1 = b(1).YEndPoints;
        labels1 = string(round(b(1).YData,2));
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')

        xtips2 = b(2).XEndPoints;
        ytips2 = b(2).YEndPoints;
        labels2 = string(round(b(2).YData,2));
        text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')
    
        xtips3 = b(3).XEndPoints;
        ytips3 = b(3).YEndPoints;
        labels3 = string(round(b(3).YData,2));
        text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
        'VerticalAlignment','bottom')

% saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\p38_NFkB_Both_', 'DoseResponses', 'ValidationAccuracies', '.fig'])
% saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\p38_NFkB_Both_', 'DoseResponses', 'ValidationAccuracies','.svg'])
%% Figure 5 Dose Responses x4, p38 only, NFkB only, Both F1 Scores by class
name =  {'P3C4_Doses_Pool','CpG_Doses_Pool','LPS_Doses_Pool','TNF_Doses_Pool'};
score_names = 'Fscore';
order = [1,2,3,4,5];
colors = setcolors;

figure('Position',[391, 156,668, 365]);
tiledlayout(2,2)
for i = 1:numel(name)
nexttile
 plot_data_ktr = StatsValidationKTR.(name{i}).(score_names)(order)';
 plot_data_nfkb = StatsValidationNFkB.(name{i}).(score_names)(order)';
 plot_data_both = StatsValidationBoth.(name{i}).(score_names)(order)';
 plot_data = [plot_data_ktr,plot_data_nfkb, plot_data_both];  
 plot_data = reshape(plot_data, [],3);

    b= bar(plot_data, 0.75, 'FaceColor','flat');
    ylabel('F1-score')
    xlabel('Dose')
    xticklabels(["2/6","3/6","4/6","5/6","6/6",])
    b(1).FaceColor = colors.p38;
    b(2).FaceColor = colors.mVen;
    b(3).FaceColor = colors.mVenKTR;

        ylim([0 1])
        title([name{i},'5-fold cross-validation'])
end
% saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\MachineLearning\p38_NFkB_Both_''DoseResponses''_F1score_Class' '.fig'])
% saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\MachineLearning\p38_NFkB_Both_''DoseResponses''_F1score_Class' '.svg'])