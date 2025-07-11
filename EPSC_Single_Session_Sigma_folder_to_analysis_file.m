% Cell Align Pipeline plots Event-Aligned Timestamps for One Session
% Derived from code: time stamp cell session and compile single session

clear; close all; clc;

cell_1_name = '2023_07_13_0001_1';
cell_2_name = '2023_07_13_0001_2';
cell_names_eq = eq(cell_1_name, cell_2_name);
cell_session = cell_1_name(cell_names_eq);
cells_co_occur = 0.1; % in seconds, co-ocurrence defined across full session 
cells_corr_thresh = 0.15;  % minimum R correlation for inclusion
cells_std_thresh = 1;   % for thresholding currents from noise

% Load .mat for EPSC Data, Timestamp, and Sample Frequency
[fnEPSC, drDECMAT, ~] = uigetfile('*.xlsx',' Pick the Excel Data file'); % also defines root folder
[pathfile,align_namefile,extfile] = fileparts([drDECMAT fnEPSC]);

disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp(['Loading Cell Pair Data ' cell_session(1:end-1)])

% load Cell Channel Data
cell_1_channel = readtable([drDECMAT fnEPSC],'Sheet',cell_1_name,'VariableNamingRule','preserve');  % Cell 1 channel
cell_2_channel = readtable([drDECMAT fnEPSC],'Sheet',cell_2_name,'VariableNamingRule','preserve');  % Cell 2 channel

disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
disp('Start and Initialize Parameters for Analysis')

events_samp_start = 50;   % 50 samples = -5 ms, total 100 Hz
events_samp_end = 50;     % 50 samples = +5 ms, total 100 Hz

cell_1_std_sigma = cells_std_thresh; % sigma threshold to accept epsc (from z-score)
cell_2_std_sigma = cells_std_thresh; % sigma threshold to accept epsc (from z-score)
debug_figures = false;
plot_figures = true;
save_figures = false;
write_tables = false;

disp('Parameters Initialized')

% %% Load Excel EPSC Data - Load Raw .xlsx for Cell Peaks

cell_1_indx = cell_1_channel.('Event Num.');    % peak event index
cell_1_time = cell_1_channel.('Event Time (s)');  % peak event time in s
cell_1_base = cell_1_channel.('Baseline (pA)');  % moving window in pA
cell_1_peak = cell_1_channel.('Peak (pA)');  % peak current from 0 in pA
cell_1_amp = cell_1_channel.('Amplitude (pA)');  % peak current from baseline in pA
cell_1_rise = cell_1_channel.('Rise Time (ms)');  % 10% to 90% in ms
cell_1_halfwidth = cell_1_channel.('Half-Width (ms)');  % rise to decay in ms
cell_1_decay = cell_1_channel.('Decay % (ms)');  % 90% to 10% in ms
cell_1_AUC = cell_1_channel.('AUC (pA ms)');  % area under rise to decay in pA*ms
cell_1_AUCtime = cell_1_channel.('AUC Time (ms)');  % in ms

cell_2_indx = cell_2_channel.('Event Num.');    % peak event index
cell_2_time = cell_2_channel.('Event Time (s)');  % peak event time in s
cell_2_base = cell_2_channel.('Baseline (pA)');  % moving window in pA
cell_2_peak = cell_2_channel.('Peak (pA)');  % peak current from 0 in pA
cell_2_amp = cell_2_channel.('Amplitude (pA)');  % peak current from baseline in pA
cell_2_rise = cell_2_channel.('Rise Time (ms)');  % 10% to 90% in ms
cell_2_halfwidth = cell_2_channel.('Half-Width (ms)');  % rise to decay in ms
cell_2_decay = cell_2_channel.('Decay % (ms)');  % 90% to 10% in ms
cell_2_AUC = cell_2_channel.('AUC (pA ms)');  % area under rise to decay in pA*ms
cell_2_AUCtime = cell_2_channel.('AUC Time (ms)');  % in ms

base_line = [cell_1_base;cell_2_base];
dataend_1 = ceil(cell_1_time(end));
dataend_2 = ceil(cell_2_time(end));
dataend_end = max(dataend_1,dataend_2);

epsc_samp_freq = 10000;
epsc_time_line = 0:1/epsc_samp_freq:dataend_end;
epsc_sess_dur = epsc_time_line(end);

cell_1_event_count = numel(cell_1_time);
cell_2_event_count = numel(cell_2_time);

disp(['MATLAB Data for ' align_namefile ' Loaded'])
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')

% %% EPSC Events - Identify Cell Activity

win_samp_width = events_samp_start + events_samp_end;
delta_t = 1/epsc_samp_freq;
epsc_time_ms = 1000*(delta_t:delta_t:(dataend_end - 2*delta_t));  % change s to ms

% for peak current thresholding in pA
base_mu_mean = mean(base_line);   % current baseline mean
base_sigma_dev = std(base_line);  % current baseline sigma

% set peak current above amplitude threshold in pA
sel_thresh = abs(cells_std_thresh.*base_sigma_dev) + abs(base_mu_mean);

cell_1_index = find(abs(cell_1_peak) >= sel_thresh);
cell_1_event_amps = cell_1_peak(cell_1_index);
cell_1_event_times = cell_1_time(cell_1_index);

cell_2_index = find(abs(cell_2_peak) >= sel_thresh);
cell_2_event_amps = cell_2_peak(cell_2_index);
cell_2_event_times = cell_2_time(cell_2_index);

disp(['EPSC Session Duration: ' num2str(epsc_sess_dur) ' s']);
disp('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
   
disp('>>>>>>>>>>>>>>>>>>>>>>>>>')
disp(['Cell 1 total: ' num2str(numel(cell_1_time)) ' events']);
disp(['Cell 1 mean epsc: ' num2str(mean(cell_1_peak)) ' pA']);
disp(['Cell 1 rate: ' num2str(round(numel(cell_1_time)/((epsc_sess_dur - cell_1_time(1))),4)) ' Hz']);
disp('>>>>>>>>>>>>>>>>>>>>>>>>>')
disp(['Cell 2 total: ' num2str(numel(cell_2_time)) ' events']);
disp(['Cell 2 mean epsc: ' num2str(mean(cell_2_peak)) ' pA']);
disp(['Cell 2 rate: ' num2str(round(numel(cell_2_time)/((epsc_sess_dur - cell_2_time(1))),4)) ' Hz']);
disp('>>>>>>>>>>>>>>>>>>>>>>>>>')

if debug_figures==true
    pause
    if save_figures==false
        close all
    end
end

% %% EPSC Events - Correlate Alignments

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
    temp_peak_2 = cell_2_peak(ievent2);
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

mean_delay_1_2 = mean(cell_to_1_ts_trim,'omitnan');
mean_delay_2_1 = mean(cell_to_2_ts_trim,'omitnan');

sum_overlaps_1_2 = sum((numel(cell_delay_trim_1)+numel(cell_delay_trim_2))/2);

% Define Cross Correlation
cell_1_ts_trim = event_cell_1_ts(~isnan(event_cell_1_ts));
cell_2_ts_trim = event_cell_2_ts(~isnan(event_cell_2_ts));
[corr_ts_wide_1_2,lags_ts_wide_1_2] = funct_cowan_corr(cell_2_ts_trim,cell_1_ts_trim,cells_co_occur,10,'none');
[corr_ts_zoom_1_2,lags_ts_zoom_1_2] = funct_cowan_corr(cell_2_ts_trim,cell_1_ts_trim,cells_co_occur/10,10,'none');
[corr_ts_narr_1_2,lags_ts_narr_1_2] = funct_cowan_corr(cell_2_ts_trim,cell_1_ts_trim,cells_co_occur/100,10,'none');

lim_corr_wide = [0 ceil(10*max(corr_ts_wide_1_2)+cells_corr_thresh)]/10;
lim_corr_zoom = [0 ceil(10*max(corr_ts_zoom_1_2)+cells_corr_thresh)]/10;
lim_corr_narr = [0 ceil(10*max(corr_ts_narr_1_2)+cells_corr_thresh)]/10;

[wide_peaks, wide_index] = findpeaks(corr_ts_wide_1_2);
max_corr_wide = max(corr_ts_wide_1_2);
max_lag_wide = max(lags_ts_wide_1_2(wide_index));
max_peak_wide = max(corr_ts_wide_1_2(wide_index));
[zoom_peaks, zoom_index] = findpeaks(corr_ts_zoom_1_2);
max_corr_zoom = max(corr_ts_zoom_1_2);
max_lag_zoom = max(lags_ts_zoom_1_2(zoom_index));
max_peak_zoom = max(corr_ts_zoom_1_2(zoom_index));
[narr_peaks, narr_index] = findpeaks(corr_ts_narr_1_2);
max_corr_narr = max(corr_ts_narr_1_2);
max_lag_narr = max(lags_ts_narr_1_2(narr_index));
max_peak_narr = max(corr_ts_narr_1_2(narr_index));

index_thresh_wide = corr_ts_wide_1_2>cells_corr_thresh;
max_thresh_wide = max(corr_ts_wide_1_2(index_thresh_wide));
index_thresh_zoom = corr_ts_zoom_1_2>cells_corr_thresh;
max_thresh_zoom = max(corr_ts_zoom_1_2(index_thresh_zoom));
index_thresh_narr = corr_ts_narr_1_2>cells_corr_thresh;
max_thresh_narr = max(corr_ts_narr_1_2(index_thresh_narr));

mean_corr_wide = mean(corr_ts_wide_1_2,'omitnan');
mean_corr_zoom = mean(corr_ts_zoom_1_2,'omitnan');
mean_corr_narr = mean(corr_ts_narr_1_2,'omitnan');

mean_lag_wide = mean(lags_ts_wide_1_2(index_thresh_wide),'omitnan');
mean_lag_zoom = mean(lags_ts_zoom_1_2(index_thresh_zoom),'omitnan');
mean_lag_narr = mean(lags_ts_narr_1_2(index_thresh_narr),'omitnan');

mean_thresh_wide = mean(corr_ts_wide_1_2(index_thresh_wide),'omitnan');
mean_thresh_zoom = mean(corr_ts_zoom_1_2(index_thresh_zoom),'omitnan');
mean_thresh_narr = mean(corr_ts_narr_1_2(index_thresh_narr),'omitnan');

Cell_1_2_Corrs_fig = figure('units', 'normalized', 'outerposition', [0 0.05 1 0.5]);
sgtitle([cell_session(1:end-1) ' Cross Correlations, ' num2str(cells_corr_thresh) ' R Threshold'],'Interpreter','none')
subplot(1,3,1)
bar(lags_ts_wide_1_2,corr_ts_wide_1_2)
hold all
scatter(lags_ts_wide_1_2(index_thresh_wide),corr_ts_wide_1_2(index_thresh_wide),10,'k','filled')
xlabel([num2str(10*cells_co_occur/1) ' s Window'])
xline(0,'--k','LineWidth',1)
xline(mean_lag_wide,'k','LineWidth',2)
ylabel('Cross Correlation (R_c_o_r_r)')
yline(cells_corr_thresh,'-.k','LineWidth',1)
ylim(lim_corr_wide)
title('Wide Cooccurrence')
subplot(1,3,2)
bar(lags_ts_zoom_1_2,corr_ts_zoom_1_2)
hold all
scatter(lags_ts_zoom_1_2(index_thresh_zoom),corr_ts_zoom_1_2(index_thresh_zoom),10,'k','filled')
xlabel([num2str(10*cells_co_occur/10) ' s Window'])
xline(0,'--k','LineWidth',1)
xline(mean_lag_zoom,'k','LineWidth',2)
ylabel('Cross Correlation (R_c_o_r_r)')
yline(cells_corr_thresh,'-.k','LineWidth',1)
ylim(lim_corr_zoom)
title('Zoomed Cooccurrence')
subplot(1,3,3)
bar(lags_ts_narr_1_2,corr_ts_narr_1_2)
hold all
scatter(lags_ts_narr_1_2(index_thresh_narr),corr_ts_narr_1_2(index_thresh_narr),10,'k','filled')
xlabel([num2str(10*cells_co_occur/100) ' s Window'])
xline(0,'--k','LineWidth',1)
xline(mean_lag_narr,'k','LineWidth',2)
ylabel('Cross Correlation (R_c_o_r_r)')
yline(cells_corr_thresh,'-.k','LineWidth',1)
ylim(lim_corr_narr)
title('Narrow Cooccurrence')

if debug_figures==true
    pause
    if save_figures==false
        close all
    end
end

pause 

close all

