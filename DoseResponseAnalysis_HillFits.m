%Script to analyse p38 and NFkB dose responses using  Hill equation fits on multiple dynamic features, importantly on percent responders

% Note: mock stimulation (dose: 0 mg/ml) represented as 100x lower than lowest stimulation dose

%Startup 
OneDrivePath = getenv('OneDrive');
%% Load full dataset here
% AllDataSets can be found at https://doi.org/10.5281/zenodo.8274567

% load([OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AnalysisScripts\AllDataSets.mat'])

%% Pool the 2 biological replicates
Stimuli = {'LPS', 'TNF','P3C4', 'CpG'};

for doses = 1:6
    for i = 1:2:8
        for f = 1:numel(fieldnames(DataSets.Data{i}(doses).metrics))
            FieldNames =fieldnames(DataSets.Data{i}(doses).metrics);
            DataSets.Pooled{(i+1)/2}(doses).metrics.(FieldNames{f}) = [DataSets.Data{i}(doses).metrics.(FieldNames{f});DataSets.Data{i+1}(doses).metrics.(FieldNames{f})];
        end 
        DataSets.Pooled{(i+1)/2}(doses).graph = DataSets.Data{i}(doses).graph;
        DataSets.Pooled{(i+1)/2}(doses).graph.var_nfkb = [DataSets.Data{i}(doses).graph.var_nfkb;DataSets.Data{i+1}(doses).graph.var_nfkb];
        DataSets.Pooled{(i+1)/2}(doses).graph.var_ktr = [DataSets.Data{i}(doses).graph.var_ktr;DataSets.Data{i+1}(doses).graph.var_ktr];
        DataSets.Pooled{(i+1)/2}(doses).graph.var_nfkb_no_base_ded = [DataSets.Data{i}(doses).graph.var_nfkb_no_base_ded;DataSets.Data{i+1}(doses).graph.var_nfkb_no_base_ded];
        DataSets.Pooled{(i+1)/2}(doses).graph.var_ktr_no_base_ded = [DataSets.Data{i}(doses).graph.var_nfkb_no_base_ded;DataSets.Data{i+1}(doses).graph.var_ktr_no_base_ded];
        DataSets.Pooled{(i+1)/2}(doses).graph.celldata = [DataSets.Data{i}(doses).graph.celldata;DataSets.Data{i+1}(doses).graph.celldata];
        DataSets.Pooled{(i+1)/2}(doses).graph.sort_metric= [DataSets.Data{i}(doses).graph.sort_metric;DataSets.Data{i+1}(doses).graph.sort_metric];
        DataSets.Pooled{(i+1)/2}(doses).graph.opt_ktr.TimeTicks = [0,2,4,6];
        DataSets.Pooled{(i+1)/2}(doses).info.dose= [DataSets.Data{i}(doses).info.dose];
        %recaclulate responder fraction from pooled
        DataSets.Pooled{(i+1)/2}(doses).metrics.responders_fraction_nfkb = nnz(DataSets.Pooled{(i+1)/2}(doses).metrics.responder_status_nfkb)/numel(DataSets.Pooled{(i+1)/2}(doses).metrics.responder_status_nfkb);
       DataSets.Pooled{(i+1)/2}(doses).metrics.responders_fraction_ktr = nnz(DataSets.Pooled{(i+1)/2}(doses).metrics.responder_status_ktr)/numel(DataSets.Pooled{(i+1)/2}(doses).metrics.responder_status_ktr);
    end
end
%% Dose response curve w Hill fit, using responder fraction
ligand = {'LPS', 'TNF','P3C4','CpG'};
doses = 6;
figure('Position',[391, 156,1104, 303]);
tiledlayout(1,4)

        count = 1;
    for l = [3,4,1,2]
        dose_info = nan(1,6);
        data_for_resp_nfkb = nan(1,doses);
        data_for_resp_ktr = nan(1,doses);
        for d = 1:doses
            data_for_resp_nfkb(d) = DataSets.Pooled{l}(d).metrics.responders_fraction_nfkb;
            data_for_resp_ktr(d) = DataSets.Pooled{l}(d).metrics.responders_fraction_ktr;
            dose_info(d) = str2double(DataSets.Pooled{l}(d).info.dose{:});
        end
    
       f_nfkb = fit(dose_info', data_for_resp_nfkb', '(maxA*x.^HillC)./(halfA_conc.^HillC + x.^HillC) + intersect','Upper',[10,max(dose_info),0.9,1], 'Lower', [0.1, min(dose_info),0, 0] );
        f_ktr = fit(dose_info', data_for_resp_ktr', '(maxA*x.^HillC)./(halfA_conc.^HillC + x.^HillC) + intersect','Upper',[10,max(dose_info),0.9,1], 'Lower', [0.1, min(dose_info),0, 0] );
        x_for_fit = logspace(log10(min(dose_info)), log10(max(dose_info)), 100);
    
    %xy plot NFkB, KTR responders
        nexttile
    
        s1 = scatter(dose_info,data_for_resp_nfkb, 'o');
        s1.MarkerEdgeColor = [0.9290 0.6940 0.1250];
        s1.MarkerFaceColor = [0.9290 0.6940 0.1250];
        
        hold on
        s2 = scatter(dose_info,data_for_resp_ktr, '*');
        s2.MarkerEdgeColor = [0, 0.4470, 0.7410];
        s2.MarkerFaceColor = [0, 0.4470, 0.7410];
        
        p1= plot(x_for_fit,f_nfkb(x_for_fit), 'Color', [0.9290 0.6940 0.1250]);
        p2 = plot(x_for_fit,f_ktr(x_for_fit), 'Color', [0, 0.4470, 0.7410]);
        
        ylim([0 1])
        xlim([dose_info(1), dose_info(doses)])
        set(gca, 'XTick', dose_info,'XTickLabels', dose_info);
        set(gca, 'xscale', 'log')
        title(['Responding cells ', ligand{l}])
        ylabel('Responder Fraction')
        xlabel('Dose')
        
        equ_nfkb = ['NFkB: (' num2str(f_nfkb.maxA) '*x.^' num2str(f_nfkb.HillC) ')./(' num2str(f_nfkb.halfA_conc) '.^' num2str(f_nfkb.HillC) '+ x.^' num2str(f_nfkb.HillC) ') +' num2str(f_nfkb.intersect)];
        equ_ktr = ['KTR: (' num2str(f_ktr.maxA) '*x.^' num2str(f_ktr.HillC) ')./(' num2str(f_ktr.halfA_conc) '.^' num2str(f_ktr.HillC) '+ x.^' num2str(f_ktr.HillC) ') +' num2str(f_ktr.intersect)];
        
        hold off
         
        HillC_nfkb(count) = f_nfkb.HillC;
        halfA_conc_nfkb(count) = f_nfkb.halfA_conc;
        HillC_ktr(count) = f_ktr.HillC;
        halfA_conc_ktr(count)= f_ktr.halfA_conc;
        count = count+1;
    end

%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Microscopy\HILL_fits_ResponderFraction', '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\Microscopy\HILL_fits_ResponderFraction','.svg'])

%% Dose response curve w Hill fit, using more metrics

ligand = {'LPS', 'TNF','P3C4','CpG'};
doses = 6;
figure('Position',[391, 156,1067, 1074]);

metrics_ktr = {'max_derivative_pk1_ktr', 'min_derivative_pk1_ktr', 'pk1_amp_ktr', 'max_amp_ktr', 'duration_ktr', 'max_pos_integral_ktr'};
metrics_nfkb= { 'max_derivative_pk1_nfkb', 'min_derivative_pk1_nfkb', 'pk1_amp_nfkb', 'max_amp_nfkb', 'duration_nfkb', 'max_pos_integral_nfkb'};
index = ([1,1,1,1,2,1]);
tiledlayout(6,4)
HillC_nfkb = nan(6,4);
halfA_conc_nfkb = nan(6,4);
HillC_ktr = nan(6,4);
halfA_conc_ktr= nan(6,4);

for m = 1:numel(metrics_ktr)    

    count = 1;
    for l = [3,4,1,2]
        dose_info = nan(1,doses);
        data_for_resp_nfkb = nan(1,doses);
        data_for_resp_ktr = nan(1,doses);
        for d = 1:doses
            %Mean
            data_for_resp_nfkb(d) = mean(DataSets.Pooled{l}(d).metrics.(metrics_nfkb{m})(:, index(m)), 'omitnan');
            data_for_resp_ktr(d) = mean(DataSets.Pooled{l}(d).metrics.(metrics_ktr{m})(:, index(m)), 'omitnan');
            dose_info(d) = str2double(DataSets.Pooled{l}(d).info.dose{:});
        end


        if m == 1||m == 3||m == 4 ||m == 6
            data_for_resp_nfkb = data_for_resp_nfkb./(max(data_for_resp_nfkb));
            data_for_resp_ktr = data_for_resp_ktr./(max(data_for_resp_ktr));
        end
        if m == 2
            data_for_resp_nfkb = data_for_resp_nfkb./(min(data_for_resp_nfkb));
            data_for_resp_ktr = data_for_resp_ktr./(min(data_for_resp_ktr));
        end
        f_nfkb = fit(dose_info', data_for_resp_nfkb', '(maxA*x.^HillC)./(halfA_conc.^HillC + x.^HillC) + intersect','Upper',[10,max(dose_info),+Inf, +Inf], 'Lower', [0.1, min(dose_info),-Inf, -Inf] );
        f_ktr = fit(dose_info', data_for_resp_ktr', '(maxA*x.^HillC)./(halfA_conc.^HillC + x.^HillC) + intersect','Upper',[10,max(dose_info),+Inf,+Inf], 'Lower', [0.1, min(dose_info),-Inf, -Inf] );
        x_for_fit = logspace(log10(min(dose_info)), log10(max(dose_info)), 100);
    
    %xy plot NFkB, KTR responders
        nexttile
        
        s1 = scatter(dose_info,data_for_resp_nfkb, 'o');
        s1.MarkerEdgeColor = [0.9290 0.6940 0.1250];
        s1.MarkerFaceColor = [0.9290 0.6940 0.1250];
        
        hold on
        s2 = scatter(dose_info,data_for_resp_ktr, '*');
        s2.MarkerEdgeColor = [0, 0.4470, 0.7410];
        s2.MarkerFaceColor = [0, 0.4470, 0.7410];
        
        p1= plot(x_for_fit,f_nfkb(x_for_fit), 'Color', [0.9290 0.6940 0.1250]);
        p2 = plot(x_for_fit,f_ktr(x_for_fit), 'Color', [0, 0.4470, 0.7410]);
        
        xlim([dose_info(1), dose_info(doses)])
       
        if m == numel(metrics_ktr); set(gca, 'XTick', dose_info,'XTickLabels', dose_info); xlabel('Dose'); else set(gca, 'XTickLabels', '');end
        set(gca, 'xscale', 'log')

        if m == 1; title([ligand{l}]); end
        labels = {'max_derivative_pk1', 'min_derivative_pk1', 'pk1_amp', 'max_amp', 'duration', 'max_pos_integral'};
        if l == 3; ylabel(labels{m});end
            
        equ_nfkb = ['NFkB: (' num2str(f_nfkb.maxA) '*x.^' num2str(f_nfkb.HillC) ')./(' num2str(f_nfkb.halfA_conc) '.^' num2str(f_nfkb.HillC) '+ x.^' num2str(f_nfkb.HillC) ') +' num2str(f_nfkb.intersect)];
        equ_ktr = ['KTR: (' num2str(f_ktr.maxA) '*x.^' num2str(f_ktr.HillC) ')./(' num2str(f_ktr.halfA_conc) '.^' num2str(f_ktr.HillC) '+ x.^' num2str(f_ktr.HillC) ') +' num2str(f_ktr.intersect)];
        
        hold off
        
        HillC_nfkb(m,count) = f_nfkb.HillC;
        halfA_conc_nfkb(m,count) = f_nfkb.halfA_conc;
        HillC_ktr(m,count) = f_ktr.HillC;
        halfA_conc_ktr(m,count)= f_ktr.halfA_conc;
        count = count+1;
    end
end
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\Writing\p38 paper\AllFigure\Microscopy\HILL_fits_Metrics', '.fig'])
%saveas(gcf, [OneDrivePath, '\PostDoc UCLA\1 Post Doc UCLA\\Writing\p38 paper\AllFigure\Microscopy\HILL_fits_Metrics','.svg'])

%%
Fold_KTRvsNFkB = halfA_conc_ktr./halfA_conc_nfkb;