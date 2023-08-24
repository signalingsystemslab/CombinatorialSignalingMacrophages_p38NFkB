% Script to train and evaluate decision tree ensemble machine learning
% models for classification of NFkB activity trajectories to stimulus
% identity and dose

% Input data availabe as csv here:
% in this repository: 'CombinatorialSignalingMacrophages_p38NFkB\DecisionTreeClassifiers\Inputs_2PoolRep_DecTrees' 

%Startup 
OneDrivePath = getenv('OneDrive');
%% All 4 Stim, Dose 6
trainingCondition = 'CC_input_AllLig_Rep1_Dose6_nfkb';
name = 'AllStimHighDosePool';
[ModelNFkB.(name), validationAccuracyNFkB.(name),validationPredictionsNFkB.(name), validationScoresNFkB.(name), validationAccuracyNFkB2.(name), StatsValidationNFkB.(name)] = Tree_Ensemble_Classification(trainingCondition,name, [1:4]');
%% PAMPs only, 3 Stim, Dose 6
trainingCondition = 'CC_input_AllLig_Rep1_Dose6_nfkb';
name = 'PAMPsHighDosePool';
[ModelNFkB.(name), validationAccuracyNFkB.(name),validationPredictionsNFkB.(name), validationScoresNFkB.(name), validationAccuracyNFkB2.(name), StatsValidationNFkB.(name)] = Tree_Ensemble_Classification(trainingCondition,name, [1,3,4]');
%% Dose responses - Doses 2 to 6
DosesToClassify = [2:6];

% LPS Pool
trainingCondition = 'CC_input_LPS_Rep1_AllDoses_nfkb';
name = 'LPS_Doses_Pool';
[ModelNFkB.(name), validationAccuracyNFkB.(name),validationPredictionsNFkB.(name), validationScoresNFkB.(name), validationAccuracyNFkB2.(name), StatsValidationNFkB.(name)] = Tree_Ensemble_Classification(trainingCondition, name, DosesToClassify');

% TNF Pool
trainingCondition = 'CC_input_TNF_Rep1_AllDoses_nfkb';
name = 'TNF_Doses_Pool';
[ModelNFkB.(name), validationAccuracyNFkB.(name),validationPredictionsNFkB.(name), validationScoresNFkB.(name), validationAccuracyNFkB2.(name), StatsValidationNFkB.(name)] = Tree_Ensemble_Classification(trainingCondition, name, DosesToClassify');

% P3C4 Pool
trainingCondition = 'CC_input_P3C4_Rep1_AllDoses_nfkb';
name = 'P3C4_Doses_Pool';
[ModelNFkB.(name), validationAccuracyNFkB.(name),validationPredictionsNFkB.(name), validationScoresNFkB.(name), validationAccuracyNFkB2.(name), StatsValidationNFkB.(name)] = Tree_Ensemble_Classification(trainingCondition,name, DosesToClassify');

% CpG Pool
trainingCondition = 'CC_input_CpG_Rep1_AllDoses_nfkb';
name = 'CpG_Doses_Pool';
[ModelNFkB.(name), validationAccuracyNFkB.(name),validationPredictionsNFkB.(name), validationScoresNFkB.(name), validationAccuracyNFkB2.(name), StatsValidationNFkB.(name)] = Tree_Ensemble_Classification(trainingCondition,name, DosesToClassify');

%% Function to retrieve input data (dynamic features) z-score them, randomly select cells for equal cell numbers across classes, call training + evaluation function, generate evaluation metrics
function [Model, validationAccuracy,validationPredictions, validationScores,validationAccuracy2,Stats_Validation] = Tree_Ensemble_Classification(trainingCondition,name,ClassNames)
OneDrivePath = getenv('OneDrive');
ClassNameTable = array2table(ClassNames);
ClassNameTable.Properties.VariableNames = {'Var229'};
myzscore = @(x) (x-mean(x, 'omitnan'))./std(x, 'omitnan');

% load input data (csv file, can be found here ???)
% in this repository: 'CombinatorialSignalingMacrophages_p38NFkB\DecisionTreeClassifiers\Inputs_2PoolRep_DecTrees' 
%trainingData = readtable([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Information theory\Jetka_MI_code\DataFramesPooled\', trainingCondition,'.csv']);

trainingData = trainingData(ismember( trainingData(:,1),ClassNameTable),:);

trainingDataZscored = [trainingData(:,1),array2table(myzscore(table2array(trainingData(:,2:end))))];
trainingDataZscored.Properties.VariableNames = trainingData.Properties.VariableNames ;

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
[Model, validationAccuracy,validationPredictions, validationScores] = trainClassifier(trainingDataZscored_EqualCell , ClassNames);

figure
Val_Matrix = confusionchart(Model.ClassificationEnsemble.Y,validationPredictions  );
Val_Matrix.RowSummary = 'row-normalized';
Val_Matrix.ColumnSummary = 'column-normalized';
title([name, ', 5-fold cross-val'])

%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z Score\PooledReps\NFkBmetrics\',name,'_5fVal_EqualCell' '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Matlab analysis\Classifiers\Z Score\PooledReps\NFkBmetrics\',name,'_5fVal_EqualCell' '.svg'])

Stats_Validation = confusionmatStats(Model.ClassificationEnsemble.Y,validationPredictions  );
validationAccuracy2 = nnz(Model.ClassificationEnsemble.Y==validationPredictions)/length(validationPredictions);  
end

%% Automatically generated function for model training and evaluation, edited to include additional evaluation metrics
function [trainedClassifier, validationAccuracy,validationPredictions, validationScores] = trainClassifier(trainingData, ClassNames)
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
%       accuracy score for each model.
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
% This code processes the data into the right shape for training the model.
inputTable = trainingData;
predictorNames = {'off_times_nfkb1', 'max_pos_integral_nfkb1', 'max_amp_nfkb1', 'pk1_amp_nfkb1', 'max2minDiff_nfkb1', 'responder_status_nfkb1', 'phase_diff1_nfkb1', 'phase_ratio1_nfkb1', 'phase_diff3_nfkb1', 'phase_ratio3_nfkb1', 'pk1_width_nfkb1', 'pk1_prom_nfkb1', 'pk1_time_nfkb1', 'timeUp2halfAmp_nfkb1', 'time2Max_nfkb1', 'timeUp2halfMax_nfkb1', 'max_derivative_pk1_nfkb1', 'timeDown2halfAmp_nfkb1', 'timeDown2halfMax_nfkb1', 'min_derivative_pk1_nfkb1', 'peakfreq_nfkb1', 'oscpower_nfkb1', 'oscbandwidth_nfkb1', 'rms_nfkb1', 'peak2rms_nfkb1', 'mean_movmad_nfkb1', 'mean_movvar_nfkb1', 'time2XmaxPosInt_nfkb1', 'time2XmaxPosInt_nfkb2', 'time2XmaxPosInt_nfkb3', 'oscfrac_nfkb1', 'oscfrac_nfkb2', 'oscfrac_nfkb3', 'oscfrac_nfkb4', 'duration_nfkb1', 'duration_nfkb2', 'duration_nfkb3', 'duration_nfkb4', 'duration_nfkb5', 'envelope_nfkb1', 'envelope_nfkb2', 'envelope_nfkb3', 'envelope_nfkb4', 'envelope_nfkb5', 'intwin3_nfkb1', 'intwin3_nfkb2', 'intwin3_nfkb3', 'intwin3_nfkb4', 'intwin3_nfkb5', 'duration_perc_pk1amp_nfkb1', 'duration_perc_pk1amp_nfkb2', 'duration_perc_pk1amp_nfkb3', 'duration_perc_pk1amp_nfkb4', 'duration_perc_pk1amp_nfkb5', 'envelope_perc_pk1amp_nfkb1', 'envelope_perc_pk1amp_nfkb2', 'envelope_perc_pk1amp_nfkb3', 'envelope_perc_pk1amp_nfkb4', 'envelope_perc_pk1amp_nfkb5', 'starttime_env_perc_pk1amp_nfkb1', 'starttime_env_perc_pk1amp_nfkb2', 'starttime_env_perc_pk1amp_nfkb3', 'starttime_env_perc_pk1amp_nfkb4', 'starttime_env_perc_pk1amp_nfkb5', 'intwin1_nfkb1', 'intwin1_nfkb2', 'intwin1_nfkb3', 'intwin1_nfkb4', 'intwin1_nfkb5', 'intwin1_nfkb6', 'intwin1_nfkb7', 'intwin_p5_nfkb1', 'intwin_p5_nfkb2', 'intwin_p5_nfkb3', 'intwin_p5_nfkb4', 'intwin_p5_nfkb5', 'intwin_p5_nfkb6', 'intwin_p5_nfkb7', 'intwin_p5_nfkb8', 'intwin_p5_nfkb9', 'intwin_p5_nfkb10', 'intwin_p5_nfkb11', 'intwin_p5_nfkb12', 'intwin_p5_nfkb13', 'intwin_p5_nfkb14', 'intwin_p5_nfkb15', 'derivatives_nfkb1', 'derivatives_nfkb2', 'derivatives_nfkb3', 'derivatives_nfkb4', 'derivatives_nfkb5', 'derivatives_nfkb6', 'derivatives_nfkb7', 'derivatives_nfkb8', 'derivatives_nfkb9', 'derivatives_nfkb10', 'derivatives_nfkb11', 'derivatives_nfkb12', 'derivatives_nfkb13', 'derivatives_nfkb14', 'derivatives_nfkb15', 'derivatives_nfkb16', 'derivatives_nfkb17', 'derivatives_nfkb18', 'derivatives_nfkb19', 'derivatives_nfkb20', 'derivatives_nfkb21', 'derivatives_nfkb22', 'derivatives_nfkb23', 'derivatives_nfkb24', 'derivatives_nfkb25', 'derivatives_nfkb26', 'derivatives_nfkb27', 'derivatives_nfkb28', 'derivatives_nfkb29', 'derivatives_nfkb30', 'derivatives_nfkb31', 'derivatives_nfkb32', 'derivatives_nfkb33', 'derivatives_nfkb34', 'derivatives_nfkb35', 'derivatives_nfkb36', 'derivatives_nfkb37', 'derivatives_nfkb38', 'derivatives_nfkb39', 'derivatives_nfkb40', 'derivatives_nfkb41', 'derivatives_nfkb42', 'derivatives_nfkb43', 'derivatives_nfkb44', 'derivatives_nfkb45', 'derivatives_nfkb46', 'derivatives_nfkb47', 'derivatives_nfkb48', 'time_series_ps_nfkb1', 'time_series_ps_nfkb2', 'time_series_ps_nfkb3', 'time_series_ps_nfkb4', 'time_series_ps_nfkb5', 'time_series_ps_nfkb6', 'time_series_ps_nfkb7', 'time_series_ps_nfkb8', 'time_series_ps_nfkb9', 'time_series_ps_nfkb10', 'time_series_ps_nfkb11', 'time_series_ps_nfkb12', 'time_series_ps_nfkb13', 'time_series_ps_nfkb14', 'time_series_ps_nfkb15', 'time_series_ps_nfkb16', 'time_series_ps_nfkb17', 'time_series_ps_nfkb18', 'time_series_ps_nfkb19', 'time_series_ps_nfkb20', 'time_series_ps_nfkb21', 'time_series_ps_nfkb22', 'time_series_ps_nfkb23', 'time_series_ps_nfkb24', 'time_series_ps_nfkb25', 'time_series_ps_nfkb26', 'time_series_ps_nfkb27', 'time_series_ps_nfkb28', 'time_series_ps_nfkb29', 'time_series_ps_nfkb30', 'time_series_ps_nfkb31', 'time_series_ps_nfkb32', 'time_series_ps_nfkb33', 'time_series_ps_nfkb34', 'time_series_ps_nfkb35', 'time_series_ps_nfkb36', 'time_series_ps_nfkb37', 'time_series_ps_nfkb38', 'time_series_ps_nfkb39', 'time_series_ps_nfkb40', 'time_series_ps_nfkb41', 'time_series_ps_nfkb42', 'time_series_ps_nfkb43', 'time_series_ps_nfkb44', 'time_series_ps_nfkb45', 'time_series_ps_nfkb46', 'time_series_ps_nfkb47', 'time_series_ps_nfkb48', 'time_series_ps_nfkb49', 'time_series_ps_nfkb50', 'time_series_ps_nfkb51', 'time_series_ps_nfkb52', 'time_series_ps_nfkb53', 'time_series_ps_nfkb54', 'time_series_ps_nfkb55', 'time_series_ps_nfkb56', 'time_series_ps_nfkb57', 'time_series_ps_nfkb58', 'time_series_ps_nfkb59', 'time_series_ps_nfkb60', 'time_series_ps_nfkb61', 'time_series_ps_nfkb62', 'time_series_ps_nfkb63', 'time_series_ps_nfkb64', 'time_series_ps_nfkb65', 'time_series_ps_nfkb66', 'time_series_ps_nfkb67', 'time_series_ps_nfkb68', 'time_series_ps_nfkb69', 'time_series_ps_nfkb70', 'time_series_ps_nfkb71', 'time_series_ps_nfkb72', 'time_series_ps_nfkb73', 'time_series_ps_nfkb74', 'time_series_ps_nfkb75', 'time_series_ps_nfkb76', 'time_series_ps_nfkb77', 'time_series_ps_nfkb78', 'time_series_ps_nfkb79', 'time_series_ps_nfkb80', 'time_series_ps_nfkb81', 'time_series_ps_nfkb82', 'time_series_ps_nfkb83', 'time_series_ps_nfkb84', 'time_series_ps_nfkb85', 'time_series_ps_nfkb86', 'time_series_ps_nfkb87', 'time_series_ps_nfkb88', 'time_series_ps_nfkb89', 'time_series_ps_nfkb90', 'time_series_ps_nfkb91', 'time_series_ps_nfkb92', 'time_series_ps_nfkb93', 'time_series_ps_nfkb94'};
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
trainedClassifier.RequiredVariables = {'derivatives_nfkb1', 'derivatives_nfkb10', 'derivatives_nfkb11', 'derivatives_nfkb12', 'derivatives_nfkb13', 'derivatives_nfkb14', 'derivatives_nfkb15', 'derivatives_nfkb16', 'derivatives_nfkb17', 'derivatives_nfkb18', 'derivatives_nfkb19', 'derivatives_nfkb2', 'derivatives_nfkb20', 'derivatives_nfkb21', 'derivatives_nfkb22', 'derivatives_nfkb23', 'derivatives_nfkb24', 'derivatives_nfkb25', 'derivatives_nfkb26', 'derivatives_nfkb27', 'derivatives_nfkb28', 'derivatives_nfkb29', 'derivatives_nfkb3', 'derivatives_nfkb30', 'derivatives_nfkb31', 'derivatives_nfkb32', 'derivatives_nfkb33', 'derivatives_nfkb34', 'derivatives_nfkb35', 'derivatives_nfkb36', 'derivatives_nfkb37', 'derivatives_nfkb38', 'derivatives_nfkb39', 'derivatives_nfkb4', 'derivatives_nfkb40', 'derivatives_nfkb41', 'derivatives_nfkb42', 'derivatives_nfkb43', 'derivatives_nfkb44', 'derivatives_nfkb45', 'derivatives_nfkb46', 'derivatives_nfkb47', 'derivatives_nfkb48', 'derivatives_nfkb5', 'derivatives_nfkb6', 'derivatives_nfkb7', 'derivatives_nfkb8', 'derivatives_nfkb9', 'duration_nfkb1', 'duration_nfkb2', 'duration_nfkb3', 'duration_nfkb4', 'duration_nfkb5', 'duration_perc_pk1amp_nfkb1', 'duration_perc_pk1amp_nfkb2', 'duration_perc_pk1amp_nfkb3', 'duration_perc_pk1amp_nfkb4', 'duration_perc_pk1amp_nfkb5', 'envelope_nfkb1', 'envelope_nfkb2', 'envelope_nfkb3', 'envelope_nfkb4', 'envelope_nfkb5', 'envelope_perc_pk1amp_nfkb1', 'envelope_perc_pk1amp_nfkb2', 'envelope_perc_pk1amp_nfkb3', 'envelope_perc_pk1amp_nfkb4', 'envelope_perc_pk1amp_nfkb5', 'intwin1_nfkb1', 'intwin1_nfkb2', 'intwin1_nfkb3', 'intwin1_nfkb4', 'intwin1_nfkb5', 'intwin1_nfkb6', 'intwin1_nfkb7', 'intwin3_nfkb1', 'intwin3_nfkb2', 'intwin3_nfkb3', 'intwin3_nfkb4', 'intwin3_nfkb5', 'intwin_p5_nfkb1', 'intwin_p5_nfkb10', 'intwin_p5_nfkb11', 'intwin_p5_nfkb12', 'intwin_p5_nfkb13', 'intwin_p5_nfkb14', 'intwin_p5_nfkb15', 'intwin_p5_nfkb2', 'intwin_p5_nfkb3', 'intwin_p5_nfkb4', 'intwin_p5_nfkb5', 'intwin_p5_nfkb6', 'intwin_p5_nfkb7', 'intwin_p5_nfkb8', 'intwin_p5_nfkb9', 'max2minDiff_nfkb1', 'max_amp_nfkb1', 'max_derivative_pk1_nfkb1', 'max_pos_integral_nfkb1', 'mean_movmad_nfkb1', 'mean_movvar_nfkb1', 'min_derivative_pk1_nfkb1', 'off_times_nfkb1', 'oscbandwidth_nfkb1', 'oscfrac_nfkb1', 'oscfrac_nfkb2', 'oscfrac_nfkb3', 'oscfrac_nfkb4', 'oscpower_nfkb1', 'peak2rms_nfkb1', 'peakfreq_nfkb1', 'phase_diff1_nfkb1', 'phase_diff3_nfkb1', 'phase_ratio1_nfkb1', 'phase_ratio3_nfkb1', 'pk1_amp_nfkb1', 'pk1_prom_nfkb1', 'pk1_time_nfkb1', 'pk1_width_nfkb1', 'responder_status_nfkb1', 'rms_nfkb1', 'starttime_env_perc_pk1amp_nfkb1', 'starttime_env_perc_pk1amp_nfkb2', 'starttime_env_perc_pk1amp_nfkb3', 'starttime_env_perc_pk1amp_nfkb4', 'starttime_env_perc_pk1amp_nfkb5', 'time2Max_nfkb1', 'time2XmaxPosInt_nfkb1', 'time2XmaxPosInt_nfkb2', 'time2XmaxPosInt_nfkb3', 'timeDown2halfAmp_nfkb1', 'timeDown2halfMax_nfkb1', 'timeUp2halfAmp_nfkb1', 'timeUp2halfMax_nfkb1', 'time_series_ps_nfkb1', 'time_series_ps_nfkb10', 'time_series_ps_nfkb11', 'time_series_ps_nfkb12', 'time_series_ps_nfkb13', 'time_series_ps_nfkb14', 'time_series_ps_nfkb15', 'time_series_ps_nfkb16', 'time_series_ps_nfkb17', 'time_series_ps_nfkb18', 'time_series_ps_nfkb19', 'time_series_ps_nfkb2', 'time_series_ps_nfkb20', 'time_series_ps_nfkb21', 'time_series_ps_nfkb22', 'time_series_ps_nfkb23', 'time_series_ps_nfkb24', 'time_series_ps_nfkb25', 'time_series_ps_nfkb26', 'time_series_ps_nfkb27', 'time_series_ps_nfkb28', 'time_series_ps_nfkb29', 'time_series_ps_nfkb3', 'time_series_ps_nfkb30', 'time_series_ps_nfkb31', 'time_series_ps_nfkb32', 'time_series_ps_nfkb33', 'time_series_ps_nfkb34', 'time_series_ps_nfkb35', 'time_series_ps_nfkb36', 'time_series_ps_nfkb37', 'time_series_ps_nfkb38', 'time_series_ps_nfkb39', 'time_series_ps_nfkb4', 'time_series_ps_nfkb40', 'time_series_ps_nfkb41', 'time_series_ps_nfkb42', 'time_series_ps_nfkb43', 'time_series_ps_nfkb44', 'time_series_ps_nfkb45', 'time_series_ps_nfkb46', 'time_series_ps_nfkb47', 'time_series_ps_nfkb48', 'time_series_ps_nfkb49', 'time_series_ps_nfkb5', 'time_series_ps_nfkb50', 'time_series_ps_nfkb51', 'time_series_ps_nfkb52', 'time_series_ps_nfkb53', 'time_series_ps_nfkb54', 'time_series_ps_nfkb55', 'time_series_ps_nfkb56', 'time_series_ps_nfkb57', 'time_series_ps_nfkb58', 'time_series_ps_nfkb59', 'time_series_ps_nfkb6', 'time_series_ps_nfkb60', 'time_series_ps_nfkb61', 'time_series_ps_nfkb62', 'time_series_ps_nfkb63', 'time_series_ps_nfkb64', 'time_series_ps_nfkb65', 'time_series_ps_nfkb66', 'time_series_ps_nfkb67', 'time_series_ps_nfkb68', 'time_series_ps_nfkb69', 'time_series_ps_nfkb7', 'time_series_ps_nfkb70', 'time_series_ps_nfkb71', 'time_series_ps_nfkb72', 'time_series_ps_nfkb73', 'time_series_ps_nfkb74', 'time_series_ps_nfkb75', 'time_series_ps_nfkb76', 'time_series_ps_nfkb77', 'time_series_ps_nfkb78', 'time_series_ps_nfkb79', 'time_series_ps_nfkb8', 'time_series_ps_nfkb80', 'time_series_ps_nfkb81', 'time_series_ps_nfkb82', 'time_series_ps_nfkb83', 'time_series_ps_nfkb84', 'time_series_ps_nfkb85', 'time_series_ps_nfkb86', 'time_series_ps_nfkb87', 'time_series_ps_nfkb88', 'time_series_ps_nfkb89', 'time_series_ps_nfkb9', 'time_series_ps_nfkb90', 'time_series_ps_nfkb91', 'time_series_ps_nfkb92', 'time_series_ps_nfkb93', 'time_series_ps_nfkb94'};
trainedClassifier.ClassificationEnsemble = classificationEnsemble;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2022a.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new table, T, use: \n  yfit = c.predictFcn(T) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nThe table, T, must contain the variables returned by: \n  c.RequiredVariables \nVariable formats (e.g. matrix/vector, datatype) must match the original training data. \nAdditional variables are ignored. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');

% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationEnsemble, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');
end