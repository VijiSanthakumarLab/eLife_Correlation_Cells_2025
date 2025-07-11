% Cell Align Pipeline plots Event-Aligned Timestamps for Colabeled Cells
% Derived from code: time stamp cell colabels and compile all cell colabel sessions

clear; close all; clc;

cells_co_occur = 0.1; % in seconds, co-ocurrence defined across full session
cells_corr_thresh = 0.15;  % minimum R correlation for inclusion
cells_std_thresh = 1;   % for thresholding currents from noise
epsc_samp_freq = 10000;
delta_t = 1/epsc_samp_freq;

% Load .mat for EPSC Data, Timestamp, and Sample Frequency
[fnEPSC, drDECMAT, ~] = uigetfile('*.xlsx',' Pick the Excel Data file'); % also defines root folder

disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('Initiating Cell Pair Data Sessions for Analysis')

a_list = [1,2,3,5,6,9,11,13];
b_list = [4,7,8,10,12,14,15];

alpha_list = {};
alpha_list(a_list,1) = {'ENNN'};
alpha_list(b_list,1) = {'ENEN'};

corr_names = strings;
corr_session = [];
corr_stat_path = pwd;
corr_stat_file = 'EPSC_Correlation';
corr_sheets = sheetnames([corr_stat_file '.xlsx']);

types_list = [];
types_list(1,a_list) = 1;
types_list(1,b_list) = 2;

%% Concatenate All Sessions

corr_Z_100_cat = [];
corr_Z_50_cat = [];
corr_Z_10_cat = [];

lag_Z_100_cat = [];
lag_Z_50_cat = [];
lag_Z_10_cat = [];

cell_sess_time = [];

cell_1_amps_cat = [];
cell_2_amps_cat = [];
cell_1_counts_cat = [];
cell_2_counts_cat = [];
cell_1_rates_cat = [];
cell_2_rates_cat = [];
cell_ts_group_cat = [];

% load Cell Session Data
for colabel_idx = 1:length(corr_sheets)
    curr_sess = corr_sheets(colabel_idx);
    curr_stat = readtable([corr_stat_path '\' corr_stat_file],'Sheet',curr_sess,'VariableNamingRule','preserve');
    temp_numbers = cellfun(@isnumeric,table2cell(curr_stat));
    temp_fields = fieldnames(curr_stat);
    temp_structs = {(temp_fields{(temp_numbers==0)})};
    for temp_idx = 1:length(temp_structs)
        curr_stat.(temp_structs{temp_idx}) = NaN;
    end
    corr_names(colabel_idx,:) = curr_sess;  % Session
    corr_session(colabel_idx,:) = table2array(curr_stat);  % Session

    cell_1_name = [char(corr_names(colabel_idx,:)) '_1'];
    cell_2_name = [char(corr_names(colabel_idx,:)) '_2'];
    cell_names_eq = eq(cell_1_name, cell_2_name);
    cell_session = cell_1_name(cell_names_eq);

    % load Cell Channel Data
    cell_1_channel = readtable([drDECMAT fnEPSC],'Sheet',cell_1_name,'VariableNamingRule','preserve');  % Cell 1 channel
    cell_2_channel = readtable([drDECMAT fnEPSC],'Sheet',cell_2_name,'VariableNamingRule','preserve');  % Cell 2 channel

    cell_1_time = cell_1_channel.('Event Time (s)');  % peak event time in s
    cell_1_base = cell_1_channel.('Baseline (pA)');  % moving window in pA
    cell_1_peak = cell_1_channel.('Peak (pA)');  % peak current from 0 in pA

    cell_2_time = cell_2_channel.('Event Time (s)');  % peak event time in s
    cell_2_base = cell_2_channel.('Baseline (pA)');  % moving window in pA
    cell_2_peak = cell_2_channel.('Peak (pA)');  % peak current from 0 in pA

    % Find Events by Timestamp
    delay_to_cell_1 = [];
    delay_to_cell_2 = [];
    event_cell_1_ts = [];
    event_cell_2_ts = [];

    for ievent1 = 1:length(cell_1_time)
        temp_event_1 = cell_1_time(ievent1);
        for ievent2 = 1:length(cell_2_time)
            temp_event_2 = cell_2_time(ievent2);
            temp_delay = temp_event_2 - temp_event_1;
            if -cells_co_occur<=temp_delay && temp_delay<=cells_co_occur
                delay_to_cell_1(ievent1,ievent2) = temp_delay;
                event_cell_1_ts(ievent1,ievent2) = temp_event_1;
            else
                delay_to_cell_1(ievent1,ievent2) = NaN;
                event_cell_1_ts(ievent1,ievent2) = NaN;
            end
        end
    end

    for ievent2 = 1:length(cell_2_time)
        temp_event_2 = cell_2_time(ievent2);
        for ievent1 = 1:length(cell_1_time)
            temp_event_1 = cell_1_time(ievent1);
            temp_delay = temp_event_1 - temp_event_2;
            if -cells_co_occur<=temp_delay && temp_delay<=cells_co_occur
                delay_to_cell_2(ievent2,ievent1) = temp_delay;
                event_cell_2_ts(ievent2,ievent1) = temp_event_2;
            else
                delay_to_cell_2(ievent2,ievent1) = NaN;
                event_cell_2_ts(ievent2,ievent1) = NaN;
            end
        end
    end

    cell_delay_trim_1 = delay_to_cell_1(~isnan(delay_to_cell_1));
    cell_delay_trim_2 = delay_to_cell_2(~isnan(delay_to_cell_2));
    cell_to_1_ts_trim = reshape(cell_delay_trim_1,1,numel(cell_delay_trim_1));
    cell_to_2_ts_trim = reshape(cell_delay_trim_2,1,numel(cell_delay_trim_2));
    cell_overlaps_1_2 = sum((numel(cell_delay_trim_1)+numel(cell_delay_trim_2))/2);

    % Define Cross Correlation
    cell_1_ts_trim = event_cell_1_ts(~isnan(event_cell_1_ts));
    cell_2_ts_trim = event_cell_2_ts(~isnan(event_cell_2_ts));
    [corr_ts_Z_100_1_2,lags_ts_Z_100_1_2] = funct_cowan_corr(cell_2_ts_trim,cell_1_ts_trim,cells_co_occur/10,10,'none');
    [corr_ts_Z_50_1_2,lags_ts_Z_50_1_2] = funct_cowan_corr(cell_2_ts_trim,cell_1_ts_trim,cells_co_occur/20,10,'none');
    [corr_ts_Z_10_1_2,lags_ts_Z_10_1_2] = funct_cowan_corr(cell_2_ts_trim,cell_1_ts_trim,cells_co_occur/100,10,'none');

    corr_Z_100_cat = cat(1,corr_Z_100_cat,corr_ts_Z_100_1_2);
    corr_Z_50_cat = cat(1,corr_Z_50_cat,corr_ts_Z_50_1_2);
    corr_Z_10_cat = cat(1,corr_Z_10_cat,corr_ts_Z_10_1_2);

    lag_Z_100_cat = cat(1,lag_Z_100_cat,lags_ts_Z_100_1_2);
    lag_Z_50_cat = cat(1,lag_Z_50_cat,lags_ts_Z_50_1_2);
    lag_Z_10_cat = cat(1,lag_Z_10_cat,lags_ts_Z_10_1_2);

    cell_1_event_peaks = mean(cell_1_peak,'omitnan');
    cell_2_event_peaks = mean(cell_2_peak,'omitnan');
    cell_1_event_count = numel(cell_1_base);
    cell_2_event_count = numel(cell_2_base);
    
    temp_sess_time = round(max([cell_1_time; cell_2_time]));
    cell_sess_time(colabel_idx) = temp_sess_time;
    cell_1_event_rate = cell_1_event_count/temp_sess_time;
    cell_2_event_rate = cell_2_event_count/temp_sess_time;
    
    cell_1_counts_cat = cat(1,cell_1_counts_cat,cell_1_event_count);
    cell_2_counts_cat = cat(1,cell_2_counts_cat,cell_2_event_count);
    cell_1_amps_cat = cat(1,cell_1_amps_cat,cell_1_event_peaks);
    cell_2_amps_cat = cat(1,cell_2_amps_cat,cell_2_event_peaks);
    cell_1_rates_cat = cat(1,cell_1_rates_cat,cell_1_event_rate);
    cell_2_rates_cat = cat(1,cell_2_rates_cat,cell_2_event_rate);
    
    if any(colabel_idx == a_list)
        cell_ts_group_cat = cat(1,cell_ts_group_cat,repmat('a',1,numel(cell_1_ts_trim))');
    elseif any(colabel_idx == b_list)
        cell_ts_group_cat = cat(1,cell_ts_group_cat,repmat('b',1,numel(cell_2_ts_trim))');
    end

    disp(['Loading Cell Pair Data ' cell_session(1:end-1)])
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    disp(['Cell 1 total: ' num2str(cell_1_event_count) ' events']);
    disp(['Cell 1 mean epsc: ' num2str(round(cell_1_event_peaks),3) ' pA']);
    disp(['Cell 1 rate: ' num2str(round(cell_1_event_rate,4)) ' Hz']);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>')
    disp(['Cell 2 total: ' num2str(cell_2_event_count) ' events']);
    disp(['Cell 2 mean epsc: ' num2str(round(cell_2_event_peaks),3) ' pA']);
    disp(['Cell 2 rate: ' num2str(round(cell_2_event_rate,4)) ' Hz']);
    disp('>>>>>>>>>>>>>>>>>>>>>>>>>')
end

%% Correlation Distributions Across Lag Window

disp('Statistics on Correlation Distributions Across Lag Windows')
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

lags_A_100 = mean(lag_Z_100_cat,1);
lags_A_50 = mean(lag_Z_50_cat,1);
lags_A_10 = mean(lag_Z_10_cat,1);

corr_A_100_a = corr_Z_100_cat(a_list,:)';
corr_A_50_a = corr_Z_50_cat(a_list,:)';
corr_A_10_a = corr_Z_10_cat(a_list,:)';
corr_A_100_b = corr_Z_100_cat(b_list,:)';
corr_A_50_b = corr_Z_50_cat(b_list,:)';
corr_A_10_b = corr_Z_10_cat(b_list,:)';

max_A_100_a = max(corr_A_100_a,[],'omitnan');
max_A_100_b = max(corr_A_100_b,[],'omitnan');
max_A_50_a = max(corr_A_50_a,[],'omitnan');
max_A_50_b = max(corr_A_50_b,[],'omitnan');
max_A_10_a = max(corr_A_10_a,[],'omitnan');
max_A_10_b = max(corr_A_10_b,[],'omitnan');

mean_A_100_a = mean(corr_A_100_a,2,'omitnan');
mean_A_100_b = mean(corr_A_100_b,2,'omitnan');
mean_A_50_a = mean(corr_A_50_a,2,'omitnan');
mean_A_50_b = mean(corr_A_50_b,2,'omitnan');
mean_A_10_a = mean(corr_A_10_a,2,'omitnan');
mean_A_10_b = mean(corr_A_10_b,2,'omitnan');

sigma_A_100_a = std(corr_A_100_a,0,2,'omitnan');
sigma_A_100_b = std(corr_A_100_b,0,2,'omitnan');
sigma_A_50_a = std(corr_A_50_a,0,2,'omitnan');
sigma_A_50_b = std(corr_A_50_b,0,2,'omitnan');
sigma_A_10_a = std(corr_A_10_a,0,2,'omitnan');
sigma_A_10_b = std(corr_A_10_b,0,2,'omitnan');

A_100_a_std_sqrt = sqrt(size(mean_A_100_a,1));
A_50_a_std_sqrt = sqrt(size(mean_A_50_a,1));
A_10_a_std_sqrt = sqrt(size(mean_A_10_a,1));
A_100_a_std_errs = sigma_A_100_a/A_100_a_std_sqrt;
A_50_a_std_errs = sigma_A_50_a/A_50_a_std_sqrt;
A_10_a_std_errs = sigma_A_10_a/A_10_a_std_sqrt;

A_100_b_std_sqrt = sqrt(size(mean_A_100_b,1));
A_50_b_std_sqrt = sqrt(size(mean_A_50_b,1));
A_10_b_std_sqrt = sqrt(size(mean_A_10_b,1));
A_100_b_std_errs = sigma_A_100_b/A_100_b_std_sqrt;
A_50_b_std_errs = sigma_A_50_b/A_50_b_std_sqrt;
A_10_b_std_errs = sigma_A_10_b/A_10_b_std_sqrt;

corr_Z_100_a_center = corr_Z_100_cat(a_list,6:16);
corr_Z_50_a_center = corr_Z_50_cat(a_list,6:16);
corr_Z_10_a_center = corr_Z_10_cat(a_list,6:16);
corr_Z_100_b_center = corr_Z_100_cat(b_list,6:16);
corr_Z_50_b_center = corr_Z_50_cat(b_list,6:16);
corr_Z_10_b_center = corr_Z_10_cat(b_list,6:16);

% disp('Is there higher maximum correlation (perioccurrence) for colabeled versus uncolabeled sessions?')

figure
hold all
shadedErrorBar(1000*lags_A_100,mean_A_100_a,A_100_a_std_errs,'lineProps','k');
shadedErrorBar(1000*lags_A_100,mean_A_100_b,A_100_b_std_errs,'lineProps','r');
xlim([1000*lags_A_100(1) 1000*lags_A_100(end)])
title('100 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])
figure
hold all
shadedErrorBar(1000*lags_A_50,mean_A_50_a,A_50_a_std_errs,'lineProps','k');
shadedErrorBar(1000*lags_A_50,mean_A_50_b,A_50_b_std_errs,'lineProps','g');
xlim([1000*lags_A_50(1) 1000*lags_A_50(end)])
title('50 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])
figure
hold all
shadedErrorBar(1000*lags_A_10,mean_A_10_a,A_10_a_std_errs,'lineProps','k');
shadedErrorBar(1000*lags_A_10,mean_A_10_b,A_10_b_std_errs,'lineProps','b');
xlim([1000*lags_A_10(1) 1000*lags_A_10(end)])
title('10 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])

disp('There is a trend in cooccurrence for colabeled versus uncolabeled sessions.')
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

pause

close all

disp('End Correlation Distributions Across Lag Window')
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

%% Correlation Statistics Across Sessions

disp('Statistics on Normality of Correlations, Distributions, Histograms')
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

disp('Are colabeled versus uncolabeled recording session durations different?')

cells_group = [repmat('C1',[15 1]); repmat('C2',[15 1])];
types_group = repmat(alpha_list,[2 1]);

[p_AN_sess_time,~,stats_AN_sess_time] = anova1(cell_sess_time,alpha_list,'off');
figure
multcompare(stats_AN_sess_time);
title('Durations')

disp('No difference in duration for colabeled versus uncolabeled sessions.')

disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('Are colabeled versus uncolabeled recording sessions intrinsically different?')

% figure('units', 'normalized', 'outerposition', [0 0 1 1]);

[p_amps_colabels,~,stats_amps_colabels] = anovan([cell_1_amps_cat;cell_2_amps_cat],{cells_group,types_group},'model','interaction','varnames',{'Amps','Session'},'display','off');
figure
multcompare(stats_amps_colabels,'Dimension',[2 2]);
title('Amplitudes')

[p_counts_colabels,~,stats_counts_colabels] = anovan([cell_1_counts_cat;cell_2_counts_cat],{cells_group,types_group},'model','interaction','varnames',{'Counts','Session'},'display','off');
figure
multcompare(stats_counts_colabels,'Dimension',[2 2]);
title('Counts')

[p_rates_colabels,~,stats_rates_colabels] = anovan([cell_1_rates_cat;cell_2_rates_cat],{cells_group,types_group},'model','interaction','varnames',{'Rates','Session'},'display','off');
figure
multcompare(stats_rates_colabels,'Dimension',[2 2]);
title('Rates')

disp('No difference in EPSC, Count, or Rate for colabeled versus uncolabeled sessions.')
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

pause

close all

%% Signal Detection Correlation Distributions

% Define binary response variables.
bin_res = 0.001;
bin_values = -bin_res:bin_res:0.50;

% Prediction Accuracy for True = Colabeled, at Choice Lag
disp('Prediction Accuracy for True = Colabeled Sessions at Choice Lag');
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

lag_interest = 11; % Choose range between 1 to 21 lags

% Define true false histograms. 
hist_100_a_corrs_cat = histcounts(corr_A_100_a(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_50_a_corrs_cat = histcounts(corr_A_50_a(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_10_a_corrs_cat = histcounts(corr_A_10_a(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_100_b_corrs_cat = histcounts(corr_A_100_b(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_50_b_corrs_cat = histcounts(corr_A_50_b(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_10_b_corrs_cat = histcounts(corr_A_10_b(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');

% Define true false histograms. 
Trues_cumsum_100 = cumsum(flip(hist_100_b_corrs_cat));
Falses_cumsum_100 = cumsum(flip(hist_100_a_corrs_cat));
AUC_ROC_T_F_100 = trapz(Falses_cumsum_100,Trues_cumsum_100);
Choice_T_F_100 = sqrt(2).*norminv(AUC_ROC_T_F_100);
Trues_cumsum_50 = cumsum(flip(hist_50_b_corrs_cat));
Falses_cumsum_50 = cumsum(flip(hist_50_a_corrs_cat));
AUC_ROC_T_F_50 = trapz(Falses_cumsum_50,Trues_cumsum_50);
Choice_T_F_50 = sqrt(2).*norminv(AUC_ROC_T_F_50);
Trues_cumsum_10 = cumsum(flip(hist_10_b_corrs_cat));
Falses_cumsum_10 = cumsum(flip(hist_10_a_corrs_cat));
AUC_ROC_T_F_10 = trapz(Falses_cumsum_10,Trues_cumsum_10);
Choice_T_F_10 = sqrt(2).*norminv(AUC_ROC_T_F_10);

disp(['Performance AUC at Bin ' num2str(lag_interest) ' = ' num2str(AUC_ROC_T_F_100*100) ' %']);
disp(['Performance d Prime at Bin ' num2str(lag_interest) ' = ' num2str(Choice_T_F_100) ' in SDT']);
disp(['Performance AUC at Bin ' num2str(lag_interest) ' = ' num2str(AUC_ROC_T_F_50*100) ' %']);
disp(['Performance d Prime at Bin ' num2str(lag_interest) ' = ' num2str(Choice_T_F_50) ' in SDT']);
disp(['Performance AUC at Bin ' num2str(lag_interest) ' = ' num2str(AUC_ROC_T_F_10*100) ' %']);
disp(['Performance d Prime at Bin ' num2str(lag_interest) ' = ' num2str(Choice_T_F_10) ' in SDT']);
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

figure
hold all
ax_plot(1) = plot(Falses_cumsum_100,Trues_cumsum_100,'r','LineWidth',2);
ax_plot(2) = plot(0:bin_res:1,0:bin_res:1,'--k','LineWidth',1);
lgd_plot = legend('box','off');
lgd_plot.String{1} = ['AUC = ' num2str(AUC_ROC_T_F_100)];
lgd_plot.String{2} = ['d prime = ' num2str(Choice_T_F_100)];
xlabel('False positive rate') 
ylabel('True positive rate')
title(['ROC 100 ms at Lag = ' num2str(1000*lags_A_100(lag_interest)) ' ms'])

figure
hold all
ax_plot(1) = plot(Falses_cumsum_50,Trues_cumsum_50,'g','LineWidth',2);
ax_plot(2) = plot(0:bin_res:1,0:bin_res:1,'--k','LineWidth',1);
lgd_plot = legend('box','off');
lgd_plot.String{1} = ['AUC = ' num2str(AUC_ROC_T_F_50)];
lgd_plot.String{2} = ['d prime = ' num2str(Choice_T_F_50)];
xlabel('False positive rate') 
ylabel('True positive rate')
title(['ROC 50 ms at Lag = ' num2str(1000*lags_A_50(lag_interest)) ' ms'])

pause

close all

%% Load Random Jitter Session Iterations

load('EPSC_Corr_Rand_Iter_0')

%% Correlation Aligned Distributions Across Lag Window

lags_A_100 = mean(lag_Z_100_cat,1);
lags_A_50 = mean(lag_Z_50_cat,1);
lags_A_10 = mean(lag_Z_10_cat,1);

corr_A_100_a = corr_Z_100_cat(a_list,:)';
corr_A_50_a = corr_Z_50_cat(a_list,:)';
corr_A_10_a = corr_Z_10_cat(a_list,:)';
corr_A_100_b = corr_Z_100_cat(b_list,:)';
corr_A_50_b = corr_Z_50_cat(b_list,:)';
corr_A_10_b = corr_Z_10_cat(b_list,:)';

mean_A_100_a = mean(corr_A_100_a,2,'omitnan');
mean_A_50_a = mean(corr_A_50_a,2,'omitnan');
mean_A_10_a = mean(corr_A_10_a,2,'omitnan');
sigma_A_100_a = std(corr_A_100_a,0,2,'omitnan');
sigma_A_50_a = std(corr_A_50_a,0,2,'omitnan');
sigma_A_10_a = std(corr_A_10_a,0,2,'omitnan');
A_100_a_std_sqrt = sqrt(size(mean_A_100_a,1));
A_50_a_std_sqrt = sqrt(size(mean_A_50_a,1));
A_10_a_std_sqrt = sqrt(size(mean_A_10_a,1));
A_100_a_std_errs = sigma_A_100_a/A_100_a_std_sqrt;
A_50_a_std_errs = sigma_A_50_a/A_50_a_std_sqrt;
A_10_a_std_errs = sigma_A_10_a/A_10_a_std_sqrt;

mean_A_100_b = mean(corr_A_100_b,2,'omitnan');
mean_A_50_b = mean(corr_A_50_b,2,'omitnan');
mean_A_10_b = mean(corr_A_10_b,2,'omitnan');
sigma_A_100_b = std(corr_A_100_b,0,2,'omitnan');
sigma_A_50_b = std(corr_A_50_b,0,2,'omitnan');
sigma_A_10_b = std(corr_A_10_b,0,2,'omitnan');
A_100_b_std_sqrt = sqrt(size(mean_A_100_b,1));
A_50_b_std_sqrt = sqrt(size(mean_A_50_b,1));
A_10_b_std_sqrt = sqrt(size(mean_A_10_b,1));
A_100_b_std_errs = sigma_A_100_b/A_100_b_std_sqrt;
A_50_b_std_errs = sigma_A_50_b/A_50_b_std_sqrt;
A_10_b_std_errs = sigma_A_10_b/A_10_b_std_sqrt;

center_mean_A_100_b = mean_A_100_b(11)
center_mean_A_100_a = mean_A_100_a(11)
center_mean_A_50_b = mean_A_50_b(11)
center_mean_A_50_a = mean_A_50_a(11)
center_mean_A_10_b = mean_A_10_b(11)
center_mean_A_10_a = mean_A_10_a(11)

center_std_dev_A_100_b = sigma_A_100_b(11)
center_std_dev_A_100_a = sigma_A_100_a(11)
center_std_dev_A_50_b = sigma_A_50_b(11)
center_std_dev_A_50_a = sigma_A_50_a(11)
center_std_dev_A_10_b = sigma_A_10_b(11)
center_std_dev_A_10_a = sigma_A_10_a(11)

Cell_Sess_Corrs_fig = figure('units', 'normalized', 'outerposition', [0 0.02 1 0.96]);
sgtitle(['Perioccurrence by Correlation Lags Across ' num2str(size(types_list(1,:),2)) ' Sessions'])
subplot(3,3,1)
hold all
shadedErrorBar(1000*lags_A_100,mean_A_100_b,A_100_b_std_errs,'lineProps','r');
xlim([1000*lags_A_100(1) 1000*lags_A_100(end)])
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.1,'-.k')
title('100 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])
subplot(3,3,2)
hold all
shadedErrorBar(1000*lags_A_50,mean_A_50_b,A_50_b_std_errs,'lineProps','g');
xlim([1000*lags_A_50(1) 1000*lags_A_50(end)])
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.1,'-.k')
title('50 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])
subplot(3,3,3)
hold all
shadedErrorBar(1000*lags_A_10,mean_A_10_b,A_10_b_std_errs,'lineProps','b');
xlim([1000*lags_A_10(1) 1000*lags_A_10(end)])
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.1,'-.k')
title('10 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])

subplot(3,3,4)
hold all
shadedErrorBar(1000*lags_A_100,mean_A_100_a,A_100_a_std_errs,'lineProps','k');
xlim([1000*lags_A_100(1) 1000*lags_A_100(end)])
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.1,'-.k')
title('100 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])
subplot(3,3,5)
hold all
shadedErrorBar(1000*lags_A_50,mean_A_50_a,A_50_a_std_errs,'lineProps','k');
xlim([1000*lags_A_50(1) 1000*lags_A_50(end)])
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.1,'-.k')
title('50 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])
subplot(3,3,6)
hold all
shadedErrorBar(1000*lags_A_10,mean_A_10_a,A_10_a_std_errs,'lineProps','k');
xlim([1000*lags_A_10(1) 1000*lags_A_10(end)])
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.1,'-.k')
title('10 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.25])

pause

close all

%% Correlation Jittered Distributions Across Lag Window

lags_J_100 = mean(lags_Z_100_iter,1,'omitnan');
lags_J_50 = mean(lags_Z_50_iter,1,'omitnan');
lags_J_10 = mean(lags_Z_10_iter,1,'omitnan');

mean_J_100_a = mean(mean_Z_100_a_iter,1,'omitnan')
mean_J_50_a = mean(mean_Z_50_a_iter,1,'omitnan')
mean_J_10_a = mean(mean_Z_10_a_iter,1,'omitnan')
J_100_a_std_devs = std(mean_Z_100_a_iter,0,1,'omitnan')
J_50_a_std_devs = std(mean_Z_50_a_iter,0,1,'omitnan')
J_10_a_std_devs = std(mean_Z_10_a_iter,0,1,'omitnan')
J_100_a_std_sqrt = sqrt(size(mean_Z_100_a_iter,1));
J_50_a_std_sqrt = sqrt(size(mean_Z_50_a_iter,1));
J_10_a_std_sqrt = sqrt(size(mean_Z_10_a_iter,1));
J_100_a_std_errs = J_100_a_std_devs/J_100_a_std_sqrt;
J_50_a_std_errs = J_50_a_std_devs/J_50_a_std_sqrt;
J_10_a_std_errs = J_10_a_std_devs/J_10_a_std_sqrt;

mean_J_100_b = mean(mean_Z_100_b_iter,1,'omitnan')
mean_J_50_b = mean(mean_Z_50_b_iter,1,'omitnan')
mean_J_10_b = mean(mean_Z_10_b_iter,1,'omitnan')
J_100_b_std_devs = std(mean_Z_100_b_iter,0,1,'omitnan')
J_50_b_std_devs = std(mean_Z_50_b_iter,0,1,'omitnan')
J_10_b_std_devs = std(mean_Z_10_b_iter,0,1,'omitnan')
J_100_b_std_sqrt = sqrt(size(mean_Z_100_b_iter,1));
J_50_b_std_sqrt = sqrt(size(mean_Z_50_b_iter,1));
J_10_b_std_sqrt = sqrt(size(mean_Z_10_b_iter,1));
J_100_b_std_errs = J_100_b_std_devs/J_100_b_std_sqrt;
J_50_b_std_errs = J_50_b_std_devs/J_50_b_std_sqrt;
J_10_b_std_errs = J_10_b_std_devs/J_10_b_std_sqrt;

center_mean_J_100_b = mean_J_100_b(11)
center_mean_J_100_a = mean_J_100_a(11)
center_mean_J_50_b = mean_J_50_b(11)
center_mean_J_50_a = mean_J_50_a(11)
center_mean_J_10_b = mean_J_10_b(11)
center_mean_J_10_a = mean_J_10_a(11)

center_std_dev_J_100_b = J_100_b_std_devs(11)
center_std_dev_J_100_a = J_100_a_std_devs(11)
center_std_dev_J_50_b = J_50_b_std_devs(11)
center_std_dev_J_50_a = J_50_a_std_devs(11)
center_std_dev_J_10_b = J_10_b_std_devs(11)
center_std_dev_J_10_a = J_10_a_std_devs(11)

[~,p_100_a_rand_sess,~,~] = ttest2(mean_Z_100_a_iter(:,11), corr_Z_100_cat(a_list,11),'Vartype','unequal');
[~,p_100_b_rand_sess,~,~] = ttest2(mean_Z_100_b_iter(:,11), corr_Z_100_cat(b_list,11),'Vartype','unequal');
[~,p_50_a_rand_sess,~,~] = ttest2(mean_Z_50_a_iter(:,11), corr_Z_50_cat(a_list,11),'Vartype','unequal');
[~,p_50_b_rand_sess,~,~] = ttest2(mean_Z_50_b_iter(:,11), corr_Z_50_cat(b_list,11),'Vartype','unequal');
[~,p_10_a_rand_sess,~,~] = ttest2(mean_Z_10_a_iter(:,11), corr_Z_10_cat(a_list,11),'Vartype','unequal');
[~,p_10_b_rand_sess,~,~] = ttest2(mean_Z_10_b_iter(:,11), corr_Z_10_cat(b_list,11),'Vartype','unequal');

Rand_Mean_Corrs_fig = figure('units', 'normalized', 'outerposition', [0 0.02 1 0.96]);
sgtitle(['Mean Correlations, Jittered Random ' num2str(num_iterations) ' Iterations'])
subplot(3,3,1)
hold all
shadedErrorBar(1000*lags_J_100,mean_J_100_b,J_100_b_std_devs,'lineProps','r');
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.06,'-.k')
title('100 ms Lag Jittered');xlabel('Lag (ms)');xlim([-100 100]);ylabel('Correlation');ylim([0 0.25])
subplot(3,3,2)
hold all
shadedErrorBar(1000*lags_J_50,mean_J_50_b,J_50_b_std_devs,'lineProps','g');
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.03,'-.k')
title('50 ms Lag Jittered');xlabel('Lag (ms)');xlim([-50 50]);ylabel('Correlation');ylim([0 0.25])
subplot(3,3,3)
hold all
shadedErrorBar(1000*lags_J_10,mean_J_10_b,J_10_b_std_devs,'lineProps','b');
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.01,'-.k')
title('10 ms Lag Jittered');xlabel('Lag (ms)');xlim([-10 10]);ylabel('Correlation');ylim([0 0.25])

subplot(3,3,4)
hold all
shadedErrorBar(1000*lags_J_100,mean_J_100_a,J_100_a_std_devs,'lineProps','k');
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.06,'-.k')
title('100 ms Lag Jittered');xlabel('Lag (ms)');xlim([-100 100]);ylabel('Correlation');ylim([0 0.25])
subplot(3,3,5)
hold all
shadedErrorBar(1000*lags_J_50,mean_J_50_a,J_50_a_std_devs,'lineProps','k');
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.03,'-.k')
title('50 ms Lag Jittered');xlabel('Lag (ms)');xlim([-50 50]);ylabel('Correlation');ylim([0 0.25])
subplot(3,3,6)
hold all
shadedErrorBar(1000*lags_J_10,mean_J_10_a,J_10_a_std_devs,'lineProps','k');
xline(0,'-k','LineWidth',10,'Alpha',0.1);yline(0,'-.k');yline(0.01,'-.k')
title('10 ms Lag Jittered');xlabel('Lag (ms)');xlim([-10 10]);ylabel('Correlation');ylim([0 0.25])

pause

close all

%% Correlation Aligned Jitter Strength - Plots

mean_Jitter_100 = mean([mean_Z_100_a_iter;mean_Z_100_b_iter],1,'omitnan');
mean_Jitter_50 = mean([mean_Z_50_a_iter;mean_Z_50_b_iter],1,'omitnan');
J_100_std_devs = std([mean_Z_100_a_iter;mean_Z_100_b_iter],0,1,'omitnan');
J_50_std_devs = std([mean_Z_50_a_iter;mean_Z_50_b_iter],0,1,'omitnan');
J_100_std_sqrt = sqrt(size([mean_Z_100_a_iter;mean_Z_100_b_iter],1));
J_50_std_sqrt = sqrt(size([mean_Z_50_a_iter;mean_Z_50_b_iter],1));
J_100_std_errs = J_100_std_devs/J_100_std_sqrt;
J_50_std_errs = J_50_std_devs/J_50_std_sqrt;

color_jitts = [25 25 25]./255;
color_100_b = [60 120 240]./255;
color_100_a = [20 60 100]./255;
color_50_b = [120 60 240]./255;
color_50_a = [80 40 80]./255;

max_a_matrix = 2*ones([length(max_A_100_a) 1]);
max_b_matrix = 1*ones([length(max_A_100_b) 1]);

% Correlations_Jitter = figure('units', 'normalized', 'outerposition', [0 0.16 1 0.67]);
% sgtitle('Perioccurrence by Correlation Lags Across Sessions')
% subplot(1,2,1)
figure
hold all
H_100_Aligns_B = shadedErrorBar(1000*lags_J_100,mean_A_100_b,A_100_b_std_errs);
H_100_Aligns_B.mainLine.Color = color_100_b; H_100_Aligns_B.mainLine.LineWidth = 1;
H_100_Aligns_B.patch.FaceColor = color_100_b; H_100_Aligns_B.patch.FaceAlpha = 0.5;
H_100_Aligns_B.edge(1).Color = color_100_b; H_100_Aligns_B.edge(2).Color = color_100_b;
H_100_Aligns_B.edge(1).LineWidth = 1; H_100_Aligns_B.edge(2).LineWidth = 1;
H_100_Aligns_A = shadedErrorBar(1000*lags_J_100,mean_A_100_a,A_100_a_std_errs);
H_100_Aligns_A.mainLine.Color = color_100_a; H_100_Aligns_A.mainLine.LineWidth = 1;
H_100_Aligns_A.patch.FaceColor = color_100_a; H_100_Aligns_A.patch.FaceAlpha = 0.5;
H_100_Aligns_A.edge(1).Color = color_100_a; H_100_Aligns_A.edge(2).Color = color_100_a;
H_100_Aligns_A.edge(1).LineWidth = 1; H_100_Aligns_A.edge(2).LineWidth = 1;
H_100_Jitters = shadedErrorBar(1000*lags_J_100,mean_Jitter_100,J_100_std_errs);
H_100_Jitters.mainLine.Color = color_jitts; H_100_Jitters.mainLine.LineWidth = 1;
H_100_Jitters.patch.FaceColor = color_jitts; H_100_Jitters.patch.FaceAlpha = 0.2;
H_100_Jitters.edge(1).Color = color_jitts; H_100_Jitters.edge(2).Color = color_jitts;
H_100_Jitters.edge(1).LineWidth = 1; H_100_Jitters.edge(2).LineWidth = 1;
xlim([1000*lags_J_100(1) 1000*lags_J_100(end)]);xline(0,'-k','LineWidth',10,'Alpha',0.1);xline(0,'k','LineWidth',2);yline(0.1,'-.k','LineWidth',1)
title('10 ms across 200 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.15])
legend({'Colabeled' 'Uncolabeled' 'Jittered'},'Box','off','Location','northwest')

% create smaller axes in top right, and plot on it
h = axes('Parent',gcf,'Position',[0.35 0.75 0.1 0.1]);
hold all
box on
scatter(max_a_matrix,max_A_100_a,10,'filled','MarkerFaceColor',color_100_a,'MarkerEdgeColor',color_100_a,'XJitter','randn','XJitterWidth',0.4)
scatter(max_b_matrix,max_A_100_b,10,'filled','MarkerFaceColor',color_100_b,'MarkerEdgeColor',color_100_b,'XJitter','randn','XJitterWidth',0.4)
set(gca,'XLim',[0 3],'YLim',[0 0.3])
yline(0.15,'--k','LineWidth',1)

% subplot(1,2,2)
figure
hold all
H_50_Aligns_B = shadedErrorBar(1000*lags_J_50,mean_A_50_b,A_50_b_std_errs);
H_50_Aligns_B.mainLine.Color = color_50_b; H_50_Aligns_B.mainLine.LineWidth = 1;
H_50_Aligns_B.patch.FaceColor = color_50_b; H_50_Aligns_B.patch.FaceAlpha = 0.5;
H_50_Aligns_B.edge(1).Color = color_50_b; H_50_Aligns_B.edge(2).Color = color_50_b;
H_50_Aligns_B.edge(1).LineWidth = 1; H_50_Aligns_B.edge(2).LineWidth = 1;
H_50_Aligns_A = shadedErrorBar(1000*lags_J_50,mean_A_50_a,A_50_a_std_errs);
H_50_Aligns_A.mainLine.Color = color_50_a; H_50_Aligns_A.mainLine.LineWidth = 1;
H_50_Aligns_A.patch.FaceColor = color_50_a; H_50_Aligns_A.patch.FaceAlpha = 0.5;
H_50_Aligns_A.edge(1).Color = color_50_a; H_50_Aligns_A.edge(2).Color = color_50_a;
H_50_Aligns_A.edge(1).LineWidth = 1; H_50_Aligns_A.edge(2).LineWidth = 1;
H_50_Jitters = shadedErrorBar(1000*lags_J_50,mean_Jitter_50,J_50_std_errs);
H_50_Jitters.mainLine.Color = color_jitts; H_50_Jitters.mainLine.LineWidth = 1;
H_50_Jitters.patch.FaceColor = color_jitts; H_50_Jitters.patch.FaceAlpha = 0.2;
H_50_Jitters.edge(1).Color = color_jitts; H_50_Jitters.edge(2).Color = color_jitts;
H_50_Jitters.edge(1).LineWidth = 1; H_50_Jitters.edge(2).LineWidth = 1;
xlim([1000*lags_J_50(1) 1000*lags_J_50(end)]);xline(0,'-k','LineWidth',10,'Alpha',0.1);xline(0,'k','LineWidth',2);yline(0.1,'-.k','LineWidth',1)
title('5 ms across 100 ms Lag Window');xlabel('Lag (ms)');ylabel('Correlation');ylim([0 0.15])
legend({'Colabeled' 'Uncolabeled' 'Jittered'},'Box','off','Location','northwest')

% create smaller axes in top right, and plot on it
h = axes('Parent',gcf,'Position',[0.78 0.75 0.1 0.1]);
hold all
box on
scatter(max_a_matrix,max_A_50_a,10,'filled','MarkerFaceColor',color_50_a,'MarkerEdgeColor',color_50_a,'XJitter','randn','XJitterWidth',0.4)
scatter(max_b_matrix,max_A_50_b,10,'filled','MarkerFaceColor',color_50_b,'MarkerEdgeColor',color_50_b,'XJitter','randn','XJitterWidth',0.4)
set(gca,'XLim',[0 3],'YLim',[0 0.3])
yline(0.10,'--k','LineWidth',1)

pause

close all

%% Correlation Distribution Signal Detection - Center Bin

% Define binary response variables.
bin_res = 0.001;
bin_values = -bin_res:bin_res:0.50;

% Prediction Accuracy for True = Colabeled, at Choice Lag
disp('Prediction Accuracy for True = Colabeled Sessions at Choice Lag');
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

color_jitts = [25 25 25]./255;
color_100_b = [60 120 240]./255;
color_100_a = [20 60 100]./255;
color_50_b = [120 60 240]./255;
color_50_a = [80 40 80]./255;

lag_interest = 11; % Choose range between 1 to 21 lags

% Define true false histograms. 
hist_100_a_corrs_cat = histcounts(corr_A_100_a(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_50_a_corrs_cat = histcounts(corr_A_50_a(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_100_b_corrs_cat = histcounts(corr_A_100_b(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');
hist_50_b_corrs_cat = histcounts(corr_A_50_b(lag_interest,:),'BinEdges',bin_values,'Normalization','probability');

% Define true false histograms. 
Trues_cumsum_100 = cumsum(flip(hist_100_b_corrs_cat));
Falses_cumsum_100 = cumsum(flip(hist_100_a_corrs_cat));
AUC_ROC_T_F_100 = trapz(Falses_cumsum_100,Trues_cumsum_100);
Choice_T_F_100 = sqrt(2).*norminv(AUC_ROC_T_F_100);
Trues_cumsum_50 = cumsum(flip(hist_50_b_corrs_cat));
Falses_cumsum_50 = cumsum(flip(hist_50_a_corrs_cat));
AUC_ROC_T_F_50 = trapz(Falses_cumsum_50,Trues_cumsum_50);
Choice_T_F_50 = sqrt(2).*norminv(AUC_ROC_T_F_50);

disp(['Performance AUC at Choice Lag = ' num2str(AUC_ROC_T_F_100*100) ' %']);
disp(['Performance d Prime at Choice Lag = ' num2str(Choice_T_F_100) ' in SDT']);
disp(['Performance AUC at Choice Lag = ' num2str(AUC_ROC_T_F_50*100) ' %']);
disp(['Performance d Prime at Choice Lag = ' num2str(Choice_T_F_50) ' in SDT']);
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

% Plot the ROC curves.
% Corrs_Prime_Choice_fig = figure('units', 'normalized', 'outerposition', [0 0.16 1 0.67]);
% subplot(1,2,1)
figure
hold all
ax_plot(1) = plot(Falses_cumsum_100,Trues_cumsum_100,'Color',color_100_b,'LineWidth',3);
ax_plot(2) = plot(0:bin_res:1,0:bin_res:1,'--k','LineWidth',2);
legend({['AUC = ' num2str(round(AUC_ROC_T_F_100,2))],'Chance = 0.50'},'box','off','Location','northwest');
xlabel('False positive rate') 
ylabel('True positive rate')
title(['Receiver Operator Characteristic 100 ms at Lag = ' num2str(1000*lags_A_100(lag_interest)) ' ms'])
% subplot(1,2,2)
figure
hold all
ax_plot(1) = plot(Falses_cumsum_50,Trues_cumsum_50,'Color',color_50_b,'LineWidth',3);
ax_plot(2) = plot(0:bin_res:1,0:bin_res:1,'--k','LineWidth',2);
legend({['AUC = ' num2str(round(AUC_ROC_T_F_50,2))],'Chance = 0.50'},'box','off','Location','northwest');
xlabel('False positive rate') 
ylabel('True positive rate')
title(['Receiver Operator Characteristic 50 ms at Lag = ' num2str(1000*lags_A_50(lag_interest)) ' ms'])

pause 

close all
