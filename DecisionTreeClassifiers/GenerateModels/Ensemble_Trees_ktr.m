% Script to train and evaluate decision tree ensemble machine learning
% models for classification of p38 activity trajectories to stimulus
% identity and dose

% Input data availabe as csv here:
% in this repository: 'CombinatorialSignalingMacrophages_p38NFkB\DecisionTreeClassifiers\Inputs_2PoolRep_DecTrees' 

%Startup 
OneDrivePath = getenv('OneDrive');
%% All 4 Stim, Dose 6
trainingCondition = 'CC_input_AllLig_Rep1_Dose6_ktr';
name = 'AllStimHighDosePool';
[ModelKTR.(name), validationAccuracyKTR.(name),validationPredictionsKTR.(name), validationScoresKTR.(name), validationAccuracyKTR2.(name), StatsValidationKTR.(name), ValidationModel.(name), CrossValidationMargins.(name)] = Tree_Ensemble_Classification(trainingCondition,name, [1:4]');
%% PAMPs only, 3 Stim, Dose 6
trainingCondition = 'CC_input_AllLig_Rep1_Dose6_ktr';
name = 'PAMPsHighDosePool';
[ModelKTR.(name), validationAccuracyKTR.(name),validationPredictionsKTR.(name), validationScoresKTR.(name), validationAccuracyKTR2.(name), StatsValidationKTR.(name)] = Tree_Ensemble_Classification(trainingCondition,name, [1,3,4]');
%% Dose responses - Doses 2 to 6
DosesToClassify = [2:6];
% LPS Pool
trainingCondition = 'CC_input_LPS_Rep1_AllDoses_ktr';
name = 'LPS_Doses_Pool';
[ModelKTR.(name), validationAccuracyKTR.(name),validationPredictionsKTR.(name), validationScoresKTR.(name), validationAccuracyKTR2.(name), StatsValidationKTR.(name)] = Tree_Ensemble_Classification(trainingCondition, name, DosesToClassify');

% TNF Pool
trainingCondition = 'CC_input_TNF_Rep1_AllDoses_ktr';
name = 'TNF_Doses_Pool';
[ModelKTR.(name), validationAccuracyKTR.(name),validationPredictionsKTR.(name), validationScoresKTR.(name), validationAccuracyKTR2.(name), StatsValidationKTR.(name)] = Tree_Ensemble_Classification(trainingCondition, name, DosesToClassify');

% P3C4 Pool
trainingCondition = 'CC_input_P3C4_Rep1_AllDoses_ktr';
name = 'P3C4_Doses_Pool';
[ModelKTR.(name), validationAccuracyKTR.(name),validationPredictionsKTR.(name), validationScoresKTR.(name), validationAccuracyKTR2.(name), StatsValidationKTR.(name)] = Tree_Ensemble_Classification(trainingCondition,name, DosesToClassify');

% CpG Pool
trainingCondition = 'CC_input_CpG_Rep1_AllDoses_ktr';
name = 'CpG_Doses_Pool';
[ModelKTR.(name), validationAccuracyKTR.(name),validationPredictionsKTR.(name), validationScoresKTR.(name), validationAccuracyKTR2.(name), StatsValidationKTR.(name)] = Tree_Ensemble_Classification(trainingCondition,name, DosesToClassify');
%% Function to retrieve input data (dynamic features) z-score them, randomly select cells for equal cell numbers across classes, call training + evaluation function, generate evaluation metrics

function [Model, validationAccuracy,validationPredictions, validationScores, validationAccuracy2,Stats_Validation,partitionedModel,CrossValidationMargins] = Tree_Ensemble_Classification(trainingCondition,name,ClassNames)
OneDrivePath = getenv('OneDrive');
ClassNameTable = array2table(ClassNames);
ClassNameTable.Properties.VariableNames = {'Var229'};
myzscore = @(x) (x-mean(x, 'omitnan'))./std(x, 'omitnan');

% load input data (csv file, can be found here:)
% in this repository: 'CombinatorialSignalingMacrophages_p38NFkB\DecisionTreeClassifiers\Inputs_2PoolRep_DecTrees' 
%trainingData = readtable([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Information theory\Jetka_MI_code\DataFramesPooled\', trainingCondition,'.csv']);
trainingData = trainingData(ismember( trainingData(:,1),ClassNameTable),:);

trainingDataZscored = [trainingData(:,1),array2table(myzscore(table2array(trainingData(:,2:end))))];
trainingDataZscored.Properties.VariableNames = trainingData.Properties.VariableNames;

% Reduce number of cells to lowest number of cells in one class, use random cells
% Uniform random selection of cells 
    trainingDataZscored = table2array(trainingDataZscored);
    [~,sorted_index] = sort(trainingDataZscored(:,1));
    trainingDataZscored = trainingDataZscored(sorted_index,:);
    classes = unique(trainingDataZscored(:,1));
    CellNumber = min(hist(trainingDataZscored(:,1),classes));   
    
    trainingDataZscored_EqualCell = [];
    for i = 1:numel(classes)
        index(i,:) = randperm(size(trainingDataZscored((trainingDataZscored(:,1)==classes(i)),:),1),CellNumber);
        trainingDataZscored_EqualCell_part(i).data = trainingDataZscored((trainingDataZscored(:,1)==classes(i)),:);
        trainingDataZscored_EqualCell_part(i).data = trainingDataZscored_EqualCell_part(i).data(index(i,:), :);
        trainingDataZscored_EqualCell = [trainingDataZscored_EqualCell; trainingDataZscored_EqualCell_part(i).data ];
    end
    trainingDataZscored_EqualCell = array2table(trainingDataZscored_EqualCell);
    trainingDataZscored_EqualCell .Properties.VariableNames = trainingData.Properties.VariableNames;
% With equal number of cells per condition
[Model, validationAccuracy,validationPredictions, validationScores,partitionedModel,CrossValidationMargins] = trainClassifier(trainingDataZscored_EqualCell , ClassNames);

figure
Val_Matrix = confusionchart(Model.ClassificationEnsemble.Y,validationPredictions  );
Val_Matrix.RowSummary = 'row-normalized';
Val_Matrix.ColumnSummary = 'column-normalized';
title([name, ', 5-fold cross-val'])

%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z Score\PooledReps\KTR metrics\',name,'_5fVal_EqualCell' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z Score\PooledReps\KTR metrics\',name,'_5fVal_EqualCell' '.svg'])
    
Stats_Validation = confusionmatStats(Model.ClassificationEnsemble.Y,validationPredictions);
validationAccuracy2 = nnz(Model.ClassificationEnsemble.Y==validationPredictions)/length(validationPredictions);  

end

%% Automatically generated function for model training and evaluation, edited to include additional evaluation metrics

function [trainedClassifier, validationAccuracy,validationPredictions, validationScores,partitionedModel,CrossValidationMargins] = trainClassifier(trainingData, ClassNames)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% Returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: A table containing the same predictor and response
%       columns as those imported into the app.
%
%  Output:
%      trainedClassifier: A struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: A function to make predictions on new
%       data.
%
%      validationAccuracy: A double containing the accuracy as a
%       percentage. In the app, the Models pane displays this overall
%       accuracy score for each ModelKTR.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a table containing at least the same predictor columns as used
% during training. For details, enter:
%   trainedClassifier.HowToPredict

% Auto-generated by MATLAB on 17-Oct-2022 15:48:16

% Extract predictors and response
% This code processes the data into the right shape for training the ModelKTR.
inputTable = trainingData;
predictorNames = {'off_times_ktr1', 'max_pos_integral_ktr1', 'max_amp_ktr1', 'pk1_amp_ktr1', 'max2minDiff_ktr1', 'responder_status_ktr1', 'phase_diff1_ktr1', 'phase_ratio1_ktr1', 'phase_diff3_ktr1', 'phase_ratio3_ktr1', 'pk1_width_ktr1', 'pk1_prom_ktr1', 'pk1_time_ktr1', 'timeUp2halfAmp_ktr1', 'time2Max_ktr1', 'timeUp2halfMax_ktr1', 'max_derivative_pk1_ktr1', 'timeDown2halfAmp_ktr1', 'timeDown2halfMax_ktr1', 'min_derivative_pk1_ktr1', 'peakfreq_ktr1', 'oscpower_ktr1', 'oscbandwidth_ktr1', 'rms_ktr1', 'peak2rms_ktr1', 'mean_movmad_ktr1', 'mean_movvar_ktr1', 'time2XmaxPosInt_ktr1', 'time2XmaxPosInt_ktr2', 'time2XmaxPosInt_ktr3', 'oscfrac_ktr1', 'oscfrac_ktr2', 'oscfrac_ktr3', 'oscfrac_ktr4', 'duration_ktr1', 'duration_ktr2', 'duration_ktr3', 'duration_ktr4', 'duration_ktr5', 'envelope_ktr1', 'envelope_ktr2', 'envelope_ktr3', 'envelope_ktr4', 'envelope_ktr5', 'intwin3_ktr1', 'intwin3_ktr2', 'intwin3_ktr3', 'intwin3_ktr4', 'intwin3_ktr5', 'duration_perc_pk1amp_ktr1', 'duration_perc_pk1amp_ktr2', 'duration_perc_pk1amp_ktr3', 'duration_perc_pk1amp_ktr4', 'duration_perc_pk1amp_ktr5', 'envelope_perc_pk1amp_ktr1', 'envelope_perc_pk1amp_ktr2', 'envelope_perc_pk1amp_ktr3', 'envelope_perc_pk1amp_ktr4', 'envelope_perc_pk1amp_ktr5', 'starttime_env_perc_pk1amp_ktr1', 'starttime_env_perc_pk1amp_ktr2', 'starttime_env_perc_pk1amp_ktr3', 'starttime_env_perc_pk1amp_ktr4', 'starttime_env_perc_pk1amp_ktr5', 'intwin1_ktr1', 'intwin1_ktr2', 'intwin1_ktr3', 'intwin1_ktr4', 'intwin1_ktr5', 'intwin1_ktr6', 'intwin1_ktr7', 'intwin_p5_ktr1', 'intwin_p5_ktr2', 'intwin_p5_ktr3', 'intwin_p5_ktr4', 'intwin_p5_ktr5', 'intwin_p5_ktr6', 'intwin_p5_ktr7', 'intwin_p5_ktr8', 'intwin_p5_ktr9', 'intwin_p5_ktr10', 'intwin_p5_ktr11', 'intwin_p5_ktr12', 'intwin_p5_ktr13', 'intwin_p5_ktr14', 'intwin_p5_ktr15', 'derivatives_ktr1', 'derivatives_ktr2', 'derivatives_ktr3', 'derivatives_ktr4', 'derivatives_ktr5', 'derivatives_ktr6', 'derivatives_ktr7', 'derivatives_ktr8', 'derivatives_ktr9', 'derivatives_ktr10', 'derivatives_ktr11', 'derivatives_ktr12', 'derivatives_ktr13', 'derivatives_ktr14', 'derivatives_ktr15', 'derivatives_ktr16', 'derivatives_ktr17', 'derivatives_ktr18', 'derivatives_ktr19', 'derivatives_ktr20', 'derivatives_ktr21', 'derivatives_ktr22', 'derivatives_ktr23', 'derivatives_ktr24', 'derivatives_ktr25', 'derivatives_ktr26', 'derivatives_ktr27', 'derivatives_ktr28', 'derivatives_ktr29', 'derivatives_ktr30', 'derivatives_ktr31', 'derivatives_ktr32', 'derivatives_ktr33', 'derivatives_ktr34', 'derivatives_ktr35', 'derivatives_ktr36', 'derivatives_ktr37', 'derivatives_ktr38', 'derivatives_ktr39', 'derivatives_ktr40', 'derivatives_ktr41', 'derivatives_ktr42', 'derivatives_ktr43', 'derivatives_ktr44', 'derivatives_ktr45', 'derivatives_ktr46', 'derivatives_ktr47', 'derivatives_ktr48', 'time_series_ps_ktr1', 'time_series_ps_ktr2', 'time_series_ps_ktr3', 'time_series_ps_ktr4', 'time_series_ps_ktr5', 'time_series_ps_ktr6', 'time_series_ps_ktr7', 'time_series_ps_ktr8', 'time_series_ps_ktr9', 'time_series_ps_ktr10', 'time_series_ps_ktr11', 'time_series_ps_ktr12', 'time_series_ps_ktr13', 'time_series_ps_ktr14', 'time_series_ps_ktr15', 'time_series_ps_ktr16', 'time_series_ps_ktr17', 'time_series_ps_ktr18', 'time_series_ps_ktr19', 'time_series_ps_ktr20', 'time_series_ps_ktr21', 'time_series_ps_ktr22', 'time_series_ps_ktr23', 'time_series_ps_ktr24', 'time_series_ps_ktr25', 'time_series_ps_ktr26', 'time_series_ps_ktr27', 'time_series_ps_ktr28', 'time_series_ps_ktr29', 'time_series_ps_ktr30', 'time_series_ps_ktr31', 'time_series_ps_ktr32', 'time_series_ps_ktr33', 'time_series_ps_ktr34', 'time_series_ps_ktr35', 'time_series_ps_ktr36', 'time_series_ps_ktr37', 'time_series_ps_ktr38', 'time_series_ps_ktr39', 'time_series_ps_ktr40', 'time_series_ps_ktr41', 'time_series_ps_ktr42', 'time_series_ps_ktr43', 'time_series_ps_ktr44', 'time_series_ps_ktr45', 'time_series_ps_ktr46', 'time_series_ps_ktr47', 'time_series_ps_ktr48', 'time_series_ps_ktr49', 'time_series_ps_ktr50', 'time_series_ps_ktr51', 'time_series_ps_ktr52', 'time_series_ps_ktr53', 'time_series_ps_ktr54', 'time_series_ps_ktr55', 'time_series_ps_ktr56', 'time_series_ps_ktr57', 'time_series_ps_ktr58', 'time_series_ps_ktr59', 'time_series_ps_ktr60', 'time_series_ps_ktr61', 'time_series_ps_ktr62', 'time_series_ps_ktr63', 'time_series_ps_ktr64', 'time_series_ps_ktr65', 'time_series_ps_ktr66', 'time_series_ps_ktr67', 'time_series_ps_ktr68', 'time_series_ps_ktr69', 'time_series_ps_ktr70', 'time_series_ps_ktr71', 'time_series_ps_ktr72', 'time_series_ps_ktr73', 'time_series_ps_ktr74', 'time_series_ps_ktr75', 'time_series_ps_ktr76', 'time_series_ps_ktr77', 'time_series_ps_ktr78', 'time_series_ps_ktr79', 'time_series_ps_ktr80', 'time_series_ps_ktr81', 'time_series_ps_ktr82', 'time_series_ps_ktr83', 'time_series_ps_ktr84', 'time_series_ps_ktr85', 'time_series_ps_ktr86', 'time_series_ps_ktr87', 'time_series_ps_ktr88', 'time_series_ps_ktr89', 'time_series_ps_ktr90', 'time_series_ps_ktr91', 'time_series_ps_ktr92', 'time_series_ps_ktr93', 'time_series_ps_ktr94'};
predictors = inputTable(:, predictorNames);
response = inputTable.Var229;
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
template = templateTree(...
    'MaxNumSplits', 20, ...
    'NumVariablesToSample', 'all');
    classificationEnsemble = fitcensemble(...
    predictors, ...
    response, ...
    'Method', 'Bag', ...
    'NumLearningCycles', 30, ...
    'Learners', template, ...
    'ClassNames', ClassNames);
% Create the result struct with predict function
predictorExtractionFcn = @(t) t(:, predictorNames);
ensemblePredictFcn = @(x) predict(classificationEnsemble, x);
trainedClassifier.predictFcn = @(x) ensemblePredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.RequiredVariables = {'derivatives_ktr1', 'derivatives_ktr10', 'derivatives_ktr11', 'derivatives_ktr12', 'derivatives_ktr13', 'derivatives_ktr14', 'derivatives_ktr15', 'derivatives_ktr16', 'derivatives_ktr17', 'derivatives_ktr18', 'derivatives_ktr19', 'derivatives_ktr2', 'derivatives_ktr20', 'derivatives_ktr21', 'derivatives_ktr22', 'derivatives_ktr23', 'derivatives_ktr24', 'derivatives_ktr25', 'derivatives_ktr26', 'derivatives_ktr27', 'derivatives_ktr28', 'derivatives_ktr29', 'derivatives_ktr3', 'derivatives_ktr30', 'derivatives_ktr31', 'derivatives_ktr32', 'derivatives_ktr33', 'derivatives_ktr34', 'derivatives_ktr35', 'derivatives_ktr36', 'derivatives_ktr37', 'derivatives_ktr38', 'derivatives_ktr39', 'derivatives_ktr4', 'derivatives_ktr40', 'derivatives_ktr41', 'derivatives_ktr42', 'derivatives_ktr43', 'derivatives_ktr44', 'derivatives_ktr45', 'derivatives_ktr46', 'derivatives_ktr47', 'derivatives_ktr48', 'derivatives_ktr5', 'derivatives_ktr6', 'derivatives_ktr7', 'derivatives_ktr8', 'derivatives_ktr9', 'duration_ktr1', 'duration_ktr2', 'duration_ktr3', 'duration_ktr4', 'duration_ktr5', 'duration_perc_pk1amp_ktr1', 'duration_perc_pk1amp_ktr2', 'duration_perc_pk1amp_ktr3', 'duration_perc_pk1amp_ktr4', 'duration_perc_pk1amp_ktr5', 'envelope_ktr1', 'envelope_ktr2', 'envelope_ktr3', 'envelope_ktr4', 'envelope_ktr5', 'envelope_perc_pk1amp_ktr1', 'envelope_perc_pk1amp_ktr2', 'envelope_perc_pk1amp_ktr3', 'envelope_perc_pk1amp_ktr4', 'envelope_perc_pk1amp_ktr5', 'intwin1_ktr1', 'intwin1_ktr2', 'intwin1_ktr3', 'intwin1_ktr4', 'intwin1_ktr5', 'intwin1_ktr6', 'intwin1_ktr7', 'intwin3_ktr1', 'intwin3_ktr2', 'intwin3_ktr3', 'intwin3_ktr4', 'intwin3_ktr5', 'intwin_p5_ktr1', 'intwin_p5_ktr10', 'intwin_p5_ktr11', 'intwin_p5_ktr12', 'intwin_p5_ktr13', 'intwin_p5_ktr14', 'intwin_p5_ktr15', 'intwin_p5_ktr2', 'intwin_p5_ktr3', 'intwin_p5_ktr4', 'intwin_p5_ktr5', 'intwin_p5_ktr6', 'intwin_p5_ktr7', 'intwin_p5_ktr8', 'intwin_p5_ktr9', 'max2minDiff_ktr1', 'max_amp_ktr1', 'max_derivative_pk1_ktr1', 'max_pos_integral_ktr1', 'mean_movmad_ktr1', 'mean_movvar_ktr1', 'min_derivative_pk1_ktr1', 'off_times_ktr1', 'oscbandwidth_ktr1', 'oscfrac_ktr1', 'oscfrac_ktr2', 'oscfrac_ktr3', 'oscfrac_ktr4', 'oscpower_ktr1', 'peak2rms_ktr1', 'peakfreq_ktr1', 'phase_diff1_ktr1', 'phase_diff3_ktr1', 'phase_ratio1_ktr1', 'phase_ratio3_ktr1', 'pk1_amp_ktr1', 'pk1_prom_ktr1', 'pk1_time_ktr1', 'pk1_width_ktr1', 'responder_status_ktr1', 'rms_ktr1', 'starttime_env_perc_pk1amp_ktr1', 'starttime_env_perc_pk1amp_ktr2', 'starttime_env_perc_pk1amp_ktr3', 'starttime_env_perc_pk1amp_ktr4', 'starttime_env_perc_pk1amp_ktr5', 'time2Max_ktr1', 'time2XmaxPosInt_ktr1', 'time2XmaxPosInt_ktr2', 'time2XmaxPosInt_ktr3', 'timeDown2halfAmp_ktr1', 'timeDown2halfMax_ktr1', 'timeUp2halfAmp_ktr1', 'timeUp2halfMax_ktr1', 'time_series_ps_ktr1', 'time_series_ps_ktr10', 'time_series_ps_ktr11', 'time_series_ps_ktr12', 'time_series_ps_ktr13', 'time_series_ps_ktr14', 'time_series_ps_ktr15', 'time_series_ps_ktr16', 'time_series_ps_ktr17', 'time_series_ps_ktr18', 'time_series_ps_ktr19', 'time_series_ps_ktr2', 'time_series_ps_ktr20', 'time_series_ps_ktr21', 'time_series_ps_ktr22', 'time_series_ps_ktr23', 'time_series_ps_ktr24', 'time_series_ps_ktr25', 'time_series_ps_ktr26', 'time_series_ps_ktr27', 'time_series_ps_ktr28', 'time_series_ps_ktr29', 'time_series_ps_ktr3', 'time_series_ps_ktr30', 'time_series_ps_ktr31', 'time_series_ps_ktr32', 'time_series_ps_ktr33', 'time_series_ps_ktr34', 'time_series_ps_ktr35', 'time_series_ps_ktr36', 'time_series_ps_ktr37', 'time_series_ps_ktr38', 'time_series_ps_ktr39', 'time_series_ps_ktr4', 'time_series_ps_ktr40', 'time_series_ps_ktr41', 'time_series_ps_ktr42', 'time_series_ps_ktr43', 'time_series_ps_ktr44', 'time_series_ps_ktr45', 'time_series_ps_ktr46', 'time_series_ps_ktr47', 'time_series_ps_ktr48', 'time_series_ps_ktr49', 'time_series_ps_ktr5', 'time_series_ps_ktr50', 'time_series_ps_ktr51', 'time_series_ps_ktr52', 'time_series_ps_ktr53', 'time_series_ps_ktr54', 'time_series_ps_ktr55', 'time_series_ps_ktr56', 'time_series_ps_ktr57', 'time_series_ps_ktr58', 'time_series_ps_ktr59', 'time_series_ps_ktr6', 'time_series_ps_ktr60', 'time_series_ps_ktr61', 'time_series_ps_ktr62', 'time_series_ps_ktr63', 'time_series_ps_ktr64', 'time_series_ps_ktr65', 'time_series_ps_ktr66', 'time_series_ps_ktr67', 'time_series_ps_ktr68', 'time_series_ps_ktr69', 'time_series_ps_ktr7', 'time_series_ps_ktr70', 'time_series_ps_ktr71', 'time_series_ps_ktr72', 'time_series_ps_ktr73', 'time_series_ps_ktr74', 'time_series_ps_ktr75', 'time_series_ps_ktr76', 'time_series_ps_ktr77', 'time_series_ps_ktr78', 'time_series_ps_ktr79', 'time_series_ps_ktr8', 'time_series_ps_ktr80', 'time_series_ps_ktr81', 'time_series_ps_ktr82', 'time_series_ps_ktr83', 'time_series_ps_ktr84', 'time_series_ps_ktr85', 'time_series_ps_ktr86', 'time_series_ps_ktr87', 'time_series_ps_ktr88', 'time_series_ps_ktr89', 'time_series_ps_ktr9', 'time_series_ps_ktr90', 'time_series_ps_ktr91', 'time_series_ps_ktr92', 'time_series_ps_ktr93', 'time_series_ps_ktr94'};
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2022a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

% Try out KFoldMargin to estimate variation
CrossValidationMargins = kfoldMargin(partitionedModel);
end