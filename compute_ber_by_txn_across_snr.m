% =========================================================================
% Compute error probability conditioned on Tx index n across all SNR files
% and plot channel-magnitude trace + temporal correlation for one SNR file.
%
% Output for Tx n = 0..3:
%   1) P(Tx = n) across all SNR files
%   2) Error probability given Tx=n
%
% Note:
%   retransInfo in this project does not include bit-error counts.
%   So this script uses "errorState" as packet error indicator.
% =========================================================================

clear; clc; close all

%% 1) Locate latest IR folder
% Set this to a specific run folder to override auto-select.
% You can use:
%   ''                                      -> auto latest
%   '2026-03-04_16-42-13_CC_64QAM_1x1_3.000000e-07'
%   '/2026-03-04_16-42-13_CC_64QAM_1x1_3.000000e-07' (leading '/' also accepted)
manual_ir_folder = '/2026-03-04_17-41-58_CC_64QAM_1x1_3.000000e-07';

% Plot controls
plot_snr_db = [];         % [] -> first available SNR file; numeric -> nearest SNR
obs_tb_count = inf;       % observation duration in TB count (set finite value, e.g. 300)
corr_max_lag = 30;        % max lag (in TB index) for temporal correlation
h_plot_in_db = true;      % true: plot |h[n]| in dB, false: linear scale

ir_root_candidates = {'./results', './result'};
ir_root = '';
for i = 1:numel(ir_root_candidates)
    if exist(ir_root_candidates{i}, 'dir')
        ir_root = ir_root_candidates{i};
        break;
    end
end
if isempty(ir_root)
    error('Neither ./results nor ./result exists.');
end

if isempty(manual_ir_folder)
    ir_folder = findLatestIrFolder(ir_root);
else
    manual_ir_folder = strtrim(manual_ir_folder);
    if startsWith(manual_ir_folder, filesep)
        manual_ir_folder = manual_ir_folder(2:end); % treat as subfolder under results
    end
    if exist(manual_ir_folder, 'dir')
        ir_folder = manual_ir_folder;
    else
        ir_folder = fullfile(ir_root, manual_ir_folder);
    end
end
if ~exist(ir_folder, 'dir')
    error('Manual IR folder not found: %s', ir_folder);
end
fprintf('IR folder: %s\n', ir_folder);

summary_file = fullfile(ir_folder, 'simulation_summary.txt');
if ~exist(summary_file, 'file')
    error('Missing summary file: %s', summary_file);
end
snr_points = parseSummarySNR(summary_file);
[summary_nframes, summary_total_slots] = parseSummaryTiming(summary_file);
fprintf('SNR points in summary: %s\n', mat2str(snr_points', 4));

%% 2) Load retrans files and accumulate by Tx n
matFiles = dir(fullfile(ir_folder, 'retrans_snr_*.mat'));
if isempty(matFiles)
    error('No retrans_snr_*.mat files in %s', ir_folder);
end

file_idx = inf(numel(matFiles), 1);
for i = 1:numel(matFiles)
    tok = regexp(matFiles(i).name, 'retrans_snr_(\d+)\.mat', 'tokens', 'once');
    if ~isempty(tok)
        file_idx(i) = str2double(tok{1});
    end
end
[file_idx, ord] = sort(file_idx);
matFiles = matFiles(ord);

target_n = 0:3;
nFiles = numel(matFiles);
n_count_by_snr = zeros(nFiles, numel(target_n));
n_err_by_snr = zeros(nFiles, numel(target_n));
total_samples_by_snr = zeros(nFiles, 1);
snr_label_by_file = NaN(nFiles, 1);
h_series_by_snr = cell(nFiles, 1);
h_source_by_snr = cell(nFiles, 1);

for i = 1:numel(matFiles)
    fpath = fullfile(ir_folder, matFiles(i).name);
    d = load(fpath);
    if ~isfield(d, 'retransInfo')
        warning('Skipping %s (missing retransInfo).', matFiles(i).name);
        continue;
    end

    entries = d.retransInfo(:);
    retrans_num = [entries.retransNum]';
    err_state = [entries.errorState]';
    [h_series_i, h_source_i] = extractChannelSeries(d);
    h_series_by_snr{i} = h_series_i;
    h_source_by_snr{i} = h_source_i;

    total_samples_by_snr(i) = numel(retrans_num);

    for k = 1:numel(target_n)
        n = target_n(k);
        idx_n = (retrans_num == n);
        n_count_by_snr(i, k) = sum(idx_n);
        n_err_by_snr(i, k) = sum(err_state(idx_n) ~= 0);
    end

    if file_idx(i) >= 1 && file_idx(i) <= numel(snr_points)
        snr_label = snr_points(file_idx(i));
    else
        snr_label = NaN;
    end
    snr_label_by_file(i) = snr_label;
    fprintf('Loaded %-20s | file_idx=%d | SNR=%.2f dB | samples=%d\n', ...
            matFiles(i).name, file_idx(i), snr_label, numel(retrans_num));
end

total_samples = sum(total_samples_by_snr);
if total_samples == 0
    error('No retransmission samples were loaded.');
end

%% 3) Compute probabilities
n_count = sum(n_count_by_snr, 1).';
n_err = sum(n_err_by_snr, 1).';
p_tx_n = n_count / total_samples;                     % P(Tx=n)
err_prob_given_n = n_err ./ max(n_count, 1);         % error probability conditioned on Tx=n
p_tx_n_by_snr = NaN(nFiles, numel(target_n));
err_prob_given_n_by_snr = NaN(nFiles, numel(target_n));
for i = 1:nFiles
    if total_samples_by_snr(i) > 0
        p_tx_n_by_snr(i, :) = n_count_by_snr(i, :) / total_samples_by_snr(i);
    end
    err_prob_given_n_by_snr(i, :) = n_err_by_snr(i, :) ./ max(n_count_by_snr(i, :), 1);
end

%% 4) Print result tables
fprintf('\n================ PER-SNR SUMMARY ================\n');
for i = 1:nFiles
    if isnan(snr_label_by_file(i))
        snr_label_str = 'N/A';
    else
        snr_label_str = sprintf('%.2f dB', snr_label_by_file(i));
    end

    fprintf('\nSNR file #%d (%s), SNR=%s, Total samples=%d\n', ...
            file_idx(i), matFiles(i).name, snr_label_str, total_samples_by_snr(i));
    fprintf('Tx n\tTotalTB(n)\tErrorCount\tErrorCount/TotalTB(n)\tP(Tx=n|SNR)\tErrorProb|Tx=n,SNR\n');
    for k = 1:numel(target_n)
        if n_count_by_snr(i, k) > 0
            err_ratio_str = sprintf('%d/%d', n_err_by_snr(i, k), n_count_by_snr(i, k));
        else
            err_ratio_str = '0/0';
        end
        fprintf('%d\t%d\t\t%d\t\t%s\t\t\t%.8f\t%.8f\n', ...
                target_n(k), n_count_by_snr(i, k), n_err_by_snr(i, k), ...
                err_ratio_str, p_tx_n_by_snr(i, k), err_prob_given_n_by_snr(i, k));
    end
    p_row = p_tx_n_by_snr(i, :);
    p_row = p_row(isfinite(p_row));
    fprintf('Sum P(Tx=n|SNR), n=0..3: %.8f\n', sum(p_row));
end

fprintf('\n================ OVERALL SUMMARY (ALL SNR) ================\n');
fprintf('\n==============================================================\n');
fprintf('Tx-index statistics across all SNR points (n = 0..3)\n');
fprintf('Note: using errorState as error indicator\n');
fprintf('--------------------------------------------------------------\n');
fprintf('Tx n\tTotalTB(n)\tErrorCount\tErrorCount/TotalTB(n)\tP(Tx=n)\t\tErrorProb|Tx=n\n');
for k = 1:numel(target_n)
    if n_count(k) > 0
        err_ratio_str = sprintf('%d/%d', n_err(k), n_count(k));
    else
        err_ratio_str = '0/0';
    end
    fprintf('%d\t%d\t\t%d\t\t%s\t\t\t%.8f\t%.8f\n', ...
            target_n(k), n_count(k), n_err(k), err_ratio_str, p_tx_n(k), err_prob_given_n(k));
end
fprintf('--------------------------------------------------------------\n');
fprintf('Sum P(Tx=n), n=0..3: %.8f\n', sum(p_tx_n));
fprintf('==============================================================\n');

%% 5) Separate plots per SNR
for i = 1:nFiles
    if isnan(snr_label_by_file(i))
        fig_name = sprintf('Tx n Stats per SNR (file_idx=%d)', file_idx(i));
        snr_title = sprintf('file\\_idx=%d', file_idx(i));
        png_name = sprintf('tx_n_stats_fileidx_%d.png', file_idx(i));
    else
        fig_name = sprintf('Tx n Stats per SNR (%.2f dB)', snr_label_by_file(i));
        snr_title = sprintf('SNR = %.2f dB', snr_label_by_file(i));
        png_name = sprintf('tx_n_stats_snr_%0.2f.png', snr_label_by_file(i));
    end

    figure('Name', fig_name, 'Position', [120 120 900 420]);
    t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile(t, 1);
    bar(target_n, p_tx_n_by_snr(i, :), 0.6);
    grid on;
    xlabel('Tx n', 'FontWeight', 'bold');
    ylabel('Probability', 'FontWeight', 'bold');
    title('P(Tx = n | SNR)', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');

    nexttile(t, 2);
    bar(target_n, err_prob_given_n_by_snr(i, :), 0.6);
    grid on;
    xlabel('Tx n', 'FontWeight', 'bold');
    ylabel('Error Probability', 'FontWeight', 'bold');
    title('Error Probability | Tx n, SNR', 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');

    title(t, sprintf('Per-SNR Tx Statistics (%s)', snr_title), 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, fullfile(ir_folder, png_name));
end

%% 6) |h[n]| and temporal correlation plot for one selected SNR
if isempty(plot_snr_db)
    selIdx = find(isfinite(snr_label_by_file), 1, 'first');
    if isempty(selIdx)
        selIdx = 1;
    end
else
    validSnrIdx = find(isfinite(snr_label_by_file));
    if isempty(validSnrIdx)
        selIdx = 1;
    else
        [~, m] = min(abs(snr_label_by_file(validSnrIdx) - plot_snr_db));
        selIdx = validSnrIdx(m);
    end
end

h_series = h_series_by_snr{selIdx};
h_source = h_source_by_snr{selIdx};
snr_sel = snr_label_by_file(selIdx);
if isnan(snr_sel)
    snr_sel_str = sprintf('file_idx=%d', file_idx(selIdx));
else
    snr_sel_str = sprintf('SNR=%.2f dB', snr_sel);
end

if isempty(h_series)
    warning('No channel series available for %s, skipping |h[n]| and correlation plots.', snr_sel_str);
else
    nObs = min(numel(h_series), max(1, floor(obs_tb_count)));
    h_obs = h_series(1:nObs);

    if isfinite(summary_nframes) && isfinite(summary_total_slots) && summary_nframes > 0 && summary_total_slots > 0
        dt_ms = 10 * summary_nframes / summary_total_slots;
        x_obs = (0:nObs-1).' * dt_ms;
        x_label = 'Time (ms)';
    else
        x_obs = (0:nObs-1).';
        x_label = 'TB index n';
    end

    h_mag = abs(h_obs);
    if h_plot_in_db
        y_mag = 20 * log10(max(h_mag, eps));
        y_label = '|h[n]| (dB)';
    else
        y_mag = h_mag;
        y_label = '|h[n]|';
    end

    figure('Name', sprintf('|h[n]| trace (%s)', snr_sel_str), 'Position', [140 120 980 420]);
    plot(x_obs, y_mag, 'LineWidth', 1.6);
    grid on;
    xlabel(x_label, 'FontWeight', 'bold');
    ylabel(y_label, 'FontWeight', 'bold');
    title(sprintf('%s | Observation TB=%d | source: %s', snr_sel_str, nObs, h_source), 'FontWeight', 'bold');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    if isnan(snr_sel)
        saveas(gcf, fullfile(ir_folder, sprintf('h_abs_trace_fileidx_%d.png', file_idx(selIdx))));
    else
        saveas(gcf, fullfile(ir_folder, sprintf('h_abs_trace_snr_%0.2f.png', snr_sel)));
    end

    [lags, rho] = temporalCorrelation(h_obs, corr_max_lag);
    if strcmp(x_label, 'Time (ms)')
        x_lag = lags * dt_ms;
        x_lag_label = 'Lag (ms)';
    else
        x_lag = lags;
        x_lag_label = 'Lag (TB)';
    end

    figure('Name', sprintf('Temporal correlation (%s)', snr_sel_str), 'Position', [160 140 980 420]);
    plot(x_lag, abs(rho), 'o-', 'LineWidth', 1.6, 'DisplayName', '|R_h(\tau)|'); hold on;
    plot(x_lag, real(rho), 's--', 'LineWidth', 1.4, 'DisplayName', 'Re\{R_h(\tau)\}');
    grid on;
    xlabel(x_lag_label, 'FontWeight', 'bold');
    ylabel('Normalized Correlation', 'FontWeight', 'bold');
    title(sprintf('%s | maxLag=%d | source: %s', snr_sel_str, min(corr_max_lag, nObs-1), h_source), 'FontWeight', 'bold');
    legend('Location', 'best');
    set(gca, 'FontSize', 12, 'FontWeight', 'bold');
    if isnan(snr_sel)
        saveas(gcf, fullfile(ir_folder, sprintf('h_temporal_corr_fileidx_%d.png', file_idx(selIdx))));
    else
        saveas(gcf, fullfile(ir_folder, sprintf('h_temporal_corr_snr_%0.2f.png', snr_sel)));
    end
end


% =========================================================================
% Helper functions
% =========================================================================

function folder = findLatestIrFolder(rootDir)
    d = dir(rootDir);
    keep = false(numel(d), 1);
    for i = 1:numel(d)
        if ~d(i).isdir
            continue;
        end
        name = d(i).name;
        if strcmp(name, '.') || strcmp(name, '..') || ~contains(name, 'IR')
            continue;
        end
        summary_ok = exist(fullfile(rootDir, name, 'simulation_summary.txt'), 'file') == 2;
        mats_ok = ~isempty(dir(fullfile(rootDir, name, 'retrans_snr_*.mat')));
        keep(i) = summary_ok && mats_ok;
    end

    d = d(keep);
    if isempty(d)
        error('No IR result folder with simulation_summary.txt and retrans_snr_*.mat under %s', rootDir);
    end

    [~, idx] = max([d.datenum]);
    folder = fullfile(rootDir, d(idx).name);
end

function snrVec = parseSummarySNR(summaryFile)
    fid = fopen(summaryFile, 'r');
    if fid < 0
        error('Cannot open file: %s', summaryFile);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = lines{1};

    snrVec = [];
    inResults = false;
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if contains(line, 'SNR (dB)')
            inResults = true;
            continue;
        end
        if ~inResults || isempty(line)
            continue;
        end

        vals = sscanf(line, '%f');
        if numel(vals) >= 3
            snrVec(end+1,1) = vals(1); %#ok<AGROW>
        end
    end

    if isempty(snrVec)
        error('No SNR values found in %s', summaryFile);
    end

    snrVec = sort(snrVec(:));
end

function [nFrames, totalSlots] = parseSummaryTiming(summaryFile)
    nFrames = NaN;
    totalSlots = NaN;

    fid = fopen(summaryFile, 'r');
    if fid < 0
        return;
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = lines{1};
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if startsWith(line, 'NFrames:')
            vals = sscanf(line, 'NFrames: %f');
            if ~isempty(vals), nFrames = vals(1); end
        elseif startsWith(line, 'Total Slots:')
            vals = sscanf(line, 'Total Slots: %f');
            if ~isempty(vals), totalSlots = vals(1); end
        end
    end
end

function [hSeries, source] = extractChannelSeries(d)
% Build instantaneous channel series for plotting.
% Only raw saved channel trace variables are accepted.

    hSeries = [];
    source = '';

    candidateVars = {'channelCoeffTrace', 'hSeries', 'hTrace'};
    for i = 1:numel(candidateVars)
        vn = candidateVars{i};
        if isfield(d, vn) && isnumeric(d.(vn)) && ~isempty(d.(vn))
            v = d.(vn)(:);
            if ~isreal(v)
                hSeries = v;
            else
                hSeries = complex(v, zeros(size(v)));
            end
            source = sprintf('raw_%s', vn);
            return;
        end
    end

    hSeries = [];
    source = 'unavailable_raw_trace_missing';
end

function [lags, rho] = temporalCorrelation(hSeries, maxLag)
% Normalized temporal autocorrelation for complex series.

    x = hSeries(:);
    x = x(isfinite(real(x)) & isfinite(imag(x)));
    if isempty(x)
        lags = 0;
        rho = complex(NaN, NaN);
        return;
    end

    x = x - mean(x);
    den = sum(abs(x).^2);
    if den <= 0
        lags = (0:maxLag).';
        rho = complex(NaN(size(lags)), NaN(size(lags)));
        return;
    end

    maxLag = min(maxLag, numel(x)-1);
    lags = (0:maxLag).';
    rho = complex(zeros(maxLag+1,1));
    for k = 0:maxLag
        rho(k+1) = sum(x(1+k:end) .* conj(x(1:end-k))) / den;
    end
end
