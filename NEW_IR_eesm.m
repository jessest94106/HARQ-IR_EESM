% =========================================================================
% EESM Link-to-System Mapping for HARQ-IR
% Uses:
%   - TDL-C HARQ-IR results from results/ (or result/)
%   - AWGN references from AWGN_results/CR1..CR4
% =========================================================================

clear; clc; close all;

%% 1. INPUT DISCOVERY
fprintf('--- Discovering Input Data ---\n');

% Optional manual override. Leave empty to auto-select newest IR run.
manual_ir_folder = './results/2026-03-04_17-41-58_CC_64QAM_1x1_3.000000e-07';

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
    ir_folder = manual_ir_folder;
end
fprintf('IR folder: %s\n', ir_folder);

awgn_root = './AWGN_results';
awgn_cr_folders = {'CR1', 'CR2', 'CR3', 'CR4'};
awgn_refs = struct('label', {}, 'snr', {}, 'bler', {});
for i = 1:numel(awgn_cr_folders)
    cr_folder = fullfile(awgn_root, awgn_cr_folders{i});
    summary_file = fullfile(cr_folder, 'simulation_summary.txt');
    if ~exist(summary_file, 'file')
        error('Missing AWGN summary: %s', summary_file);
    end
    [snr_cr, bler_cr] = parseSummarySNRBler(summary_file);
    awgn_refs(i).label = awgn_cr_folders{i};
    awgn_refs(i).snr = snr_cr;
    awgn_refs(i).bler = bler_cr;
    fprintf('AWGN %s loaded: %d points\n', awgn_cr_folders{i}, numel(snr_cr));
end

%% 2. LOAD IR RETRANSMISSION DATA
fprintf('\n--- Loading IR Retransmission Data ---\n');
ir_summary_file = fullfile(ir_folder, 'simulation_summary.txt');
if ~exist(ir_summary_file, 'file')
    error('Missing IR summary: %s', ir_summary_file);
end
[ir_snr_points, ~] = parseSummarySNRBler(ir_summary_file);
fprintf('IR SNR points loaded: %d\n', numel(ir_snr_points));

[retransInfo_ir_all, row_snr_idx] = loadRetransInfoMatrix(ir_folder);
all_retrans_num = [retransInfo_ir_all.retransNum];
max_retrans = max(all_retrans_num);
fprintf('Max retransmission index in IR data: n = %d\n', max_retrans);
fprintf('Building packet histories ...\n');
packet_hist_by_row = buildPacketHistories(retransInfo_ir_all);

if max_retrans + 1 > numel(awgn_refs)
    warning(['IR data has %d transmissions, but only %d AWGN CR folders are available. ' ...
        'Higher n values will reuse the last AWGN reference (%s).'], ...
        max_retrans + 1, numel(awgn_refs), awgn_refs(end).label);
end

%% 3. EESM FITTING LOOP
fprintf('\n--- Running EESM Mapping ---\n');

% Runtime controls
max_fit_packets = 5000;      % use subset for beta optimization only
beta_lb = 0.05;
beta_ub = 30;
use_log_bler_objective = true;
bler_floor = 1e-5;
options = optimset('MaxFunEvals',120, ...
                   'MaxIter',120, ...
                   'TolX',1e-3, ...
                   'Display','off');

results = struct('n', {}, ...
                 'awgn_ref', {}, ...
                 'beta_static', {}, ...
                 'beta_dynamic', {}, ...
                 'mse_paper', {}, ...
                 'mse_dynamic', {});
snrEff_static_by_n = cell(max_retrans + 1, 1);
snrEff_dynamic_by_n = cell(max_retrans + 1, 1);
sample_snr_idx_static_by_n = cell(max_retrans + 1, 1);
sample_snr_idx_dynamic_by_n = cell(max_retrans + 1, 1);
beta_static_by_n = NaN(max_retrans + 1, 1);
beta_dynamic_by_n = NaN(max_retrans + 1, 1);

cols = 1;
rows = 4;
subplot_width_px = 900;
subplot_height_px = 240;
figure('Name', 'HARQ-IR EESM: All Retransmissions', ...
       'Position', [100, 100, subplot_width_px*cols, subplot_height_px*rows]);

beta_static = 1.7;

for target_n = 0:max_retrans
    t_n = tic;
    awgn_idx = min(target_n + 1, numel(awgn_refs));
    SNR_awgn = awgn_refs(awgn_idx).snr;
    BLER_awgn = awgn_refs(awgn_idx).bler;
    lut = [SNR_awgn(:), BLER_awgn(:)];

    fprintf('\n--- Transmission n = %d, AWGN ref = %s ---\n', ...
        target_n, awgn_refs(awgn_idx).label);

    n_est = sum(all_retrans_num == target_n);
    sinr_cells = cell(n_est, 1);
    perStore_n = zeros(n_est, 1);
    sample_snr_idx_n = zeros(n_est, 1);
    n_keep = 0;

    for i = 1:size(retransInfo_ir_all, 1)
        row_data = retransInfo_ir_all(i, :);
        row_hist = packet_hist_by_row{i};
        for h = 1:numel(row_hist)
            retrans_nums = row_hist(h).retransNum;
            k_last = find(retrans_nums == target_n, 1, 'last');
            if isempty(k_last)
                continue;
            end
            idx = row_hist(h).idx(1:k_last);
            sinr_concat = [];
            for k = 1:numel(idx)
                sinr_re = row_data(idx(k)).sinr_per_re;
                sinr_concat = [sinr_concat, sinr_re(:)']; %#ok<AGROW>
            end

            n_keep = n_keep + 1;
            sinr_cells{n_keep} = sinr_concat;
            perStore_n(n_keep) = row_data(idx(k_last)).errorState;
            sample_snr_idx_n(n_keep) = row_snr_idx(i);
        end
    end

    if n_keep == 0
        fprintf('No IR samples for n = %d. Skipping.\n', target_n);
        continue;
    end
    sinr_cells = sinr_cells(1:n_keep);
    perStore_n = perStore_n(1:n_keep);
    sample_snr_idx_n = sample_snr_idx_n(1:n_keep);
    sinrStore_n = cellsToMatrix(sinr_cells);

    [sinrFit_n, perFit_n] = downsamplePackets(sinrStore_n, perStore_n, max_fit_packets);
    fprintf('Samples: full=%d, fit=%d\n', size(sinrStore_n,1), size(sinrFit_n,1));

    subplot(rows, cols, target_n + 1);

    if target_n == 0
        mse_func = @(beta) awgnPerSnrFittingMse(sinrFit_n, perFit_n, lut, beta, use_log_bler_objective, bler_floor);
        beta_static = fminbnd(mse_func, beta_lb, beta_ub, options);
        [binsnr, binper, snrEff_static, mse_val] = awgnPerSnrFitting(sinrStore_n, perStore_n, lut, beta_static, use_log_bler_objective, bler_floor);
        snrEff_static_by_n{target_n + 1} = snrEff_static;
        snrEff_dynamic_by_n{target_n + 1} = snrEff_static;
        valid_idx_mask_static = snrEffValidMask(sinrStore_n, perStore_n, beta_static);
        sample_snr_idx_static_by_n{target_n + 1} = sample_snr_idx_n(valid_idx_mask_static);
        sample_snr_idx_dynamic_by_n{target_n + 1} = sample_snr_idx_n(valid_idx_mask_static);
        beta_static_by_n(target_n + 1) = beta_static;
        beta_dynamic_by_n(target_n + 1) = beta_static;

        fprintf('Calibrated baseline beta = %.4f, MSE = %.6f\n', beta_static, mse_val);

        semilogy(SNR_awgn, BLER_awgn, 'k-', 'LineWidth', 2, ...
            'DisplayName', sprintf('AWGN %s', awgn_refs(awgn_idx).label)); hold on;
        semilogy(binsnr, binper, 'k*', 'MarkerSize', 7, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('Initial fit (\\beta=%.2f)', beta_static));

        results(end+1) = struct('n', target_n, ...
                                'awgn_ref', awgn_refs(awgn_idx).label, ...
                                'beta_static', beta_static, ...
                                'beta_dynamic', beta_static, ...
                                'mse_paper', mse_val, ...
                                'mse_dynamic', mse_val); %#ok<SAGROW>
    else
        [binsnr_p, binper_p, snrEff_static, mse_p] = awgnPerSnrFitting(sinrStore_n, perStore_n, lut, beta_static, use_log_bler_objective, bler_floor);

        mse_func_dyn = @(beta) awgnPerSnrFittingMse(sinrFit_n, perFit_n, lut, beta, use_log_bler_objective, bler_floor);
        beta_dyn = fminbnd(mse_func_dyn, beta_lb, beta_ub, options);
        [binsnr_d, binper_d, snrEff_dyn, mse_d] = awgnPerSnrFitting(sinrStore_n, perStore_n, lut, beta_dyn, use_log_bler_objective, bler_floor);
        snrEff_static_by_n{target_n + 1} = snrEff_static;
        snrEff_dynamic_by_n{target_n + 1} = snrEff_dyn;
        valid_idx_mask_static = snrEffValidMask(sinrStore_n, perStore_n, beta_static);
        valid_idx_mask_dynamic = snrEffValidMask(sinrStore_n, perStore_n, beta_dyn);
        sample_snr_idx_static_by_n{target_n + 1} = sample_snr_idx_n(valid_idx_mask_static);
        sample_snr_idx_dynamic_by_n{target_n + 1} = sample_snr_idx_n(valid_idx_mask_dynamic);
        beta_static_by_n(target_n + 1) = beta_static;
        beta_dynamic_by_n(target_n + 1) = beta_dyn;

        fprintf('Paper method: beta = %.4f, MSE = %.6f\n', beta_static, mse_p);
        fprintf('Dynamic beta: beta = %.4f, MSE = %.6f\n', beta_dyn, mse_d);

        semilogy(SNR_awgn, BLER_awgn, 'k-', 'LineWidth', 2, ...
            'DisplayName', sprintf('AWGN %s', awgn_refs(awgn_idx).label)); hold on;
        semilogy(binsnr_p, binper_p, 'ro', 'MarkerSize', 6, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Paper static (\\beta=%.2f)', beta_static));
        semilogy(binsnr_d, binper_d, 'b*', 'MarkerSize', 6, 'LineWidth', 1.2, ...
            'DisplayName', sprintf('Dynamic (\\beta_q=%.2f)', beta_dyn));

        results(end+1) = struct('n', target_n, ...
                                'awgn_ref', awgn_refs(awgn_idx).label, ...
                                'beta_static', beta_static, ...
                                'beta_dynamic', beta_dyn, ...
                                'mse_paper', mse_p, ...
                                'mse_dynamic', mse_d); %#ok<SAGROW>
    end

    grid on;
    ylim([1e-3 1]);
    ax = gca;
    ax.FontSize = 14;
    ax.FontWeight = 'bold';
    ax.LineWidth = 1.2;
    xlabel('Effective SNR (dB)', 'FontSize', 16, 'FontWeight', 'bold');
    ylabel('BLER', 'FontSize', 16, 'FontWeight', 'bold');
    title(sprintf('Tx n = %d (%s)', target_n, awgn_refs(awgn_idx).label), ...
          'FontSize', 17, 'FontWeight', 'bold');
    lgd = legend('Location', 'southwest');
    lgd.FontSize = 13;
    lgd.FontWeight = 'bold';
    fprintf('Elapsed for n=%d: %.2f s\n', target_n, toc(t_n));
end

saveas(gcf, 'eesm_all_retransmissions.png');

%% 4. EFFECTIVE SINR DISTRIBUTION (HISTOGRAM)
n_plot_max = min(3, max_retrans);
n_plot = 0:n_plot_max;
n_snr_points = numel(ir_snr_points);
cmap_tx = lines(numel(n_plot));
fit_static_by_snr_tx = cell(n_snr_points, numel(n_plot));
fit_dynamic_by_snr_tx = cell(n_snr_points, numel(n_plot));

vals_for_edges_static = [];
for snr_idx = 1:n_snr_points
    for n = n_plot
        vals_n = snrEff_static_by_n{n+1};
        snr_idx_n = sample_snr_idx_static_by_n{n+1};
        if isempty(vals_n) || isempty(snr_idx_n)
            continue;
        end
        valid_mask = isfinite(vals_n) & (snr_idx_n == snr_idx);
        vals_for_edges_static = [vals_for_edges_static; vals_n(valid_mask)]; %#ok<AGROW>
    end
end

if ~isempty(vals_for_edges_static)
    minEdgeStatic = floor(min(vals_for_edges_static));
    maxEdgeStatic = ceil(max(vals_for_edges_static));
    if minEdgeStatic == maxEdgeStatic
        minEdgeStatic = minEdgeStatic - 0.5;
        maxEdgeStatic = maxEdgeStatic + 0.5;
    end
    edges_static = linspace(minEdgeStatic, maxEdgeStatic, 40);

    static_cols = min(3, n_snr_points);
    static_rows = ceil(n_snr_points / static_cols);
    figure('Name', 'Effective SINR Distribution: Static Beta by SNR (Tx n=0..3)', ...
           'Position', [120, 120, 1400, 720]);
    tl_static = tiledlayout(static_rows, static_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    fprintf('\n--- Skew Generalized Normal Fits: Static Beta ---\n');
    fprintf('SNR(dB)\tTx\tMean\t\tVariance\tLambda1\t\tLambda2\n');
    fprintf('---------------------------------------------------------------------\n');
    for snr_idx = 1:n_snr_points
        ax = nexttile(tl_static, snr_idx);
        hold(ax, 'on');
        has_data = false;
        for i = 1:numel(n_plot)
            n = n_plot(i);
            vals_n = snrEff_static_by_n{n+1};
            snr_idx_n = sample_snr_idx_static_by_n{n+1};
            if isempty(vals_n) || isempty(snr_idx_n)
                continue;
            end

            valid_mask = isfinite(vals_n) & (snr_idx_n == snr_idx);
            vals = vals_n(valid_mask);
            if isempty(vals)
                continue;
            end
            fit_static = fitSkewGeneralizedNormal(vals);
            fit_static_by_snr_tx{snr_idx, i} = fit_static;
            if fit_static.ok
                fprintf('%.2f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
                        ir_snr_points(snr_idx), n, fit_static.mean, fit_static.variance, ...
                        fit_static.lambda1, fit_static.lambda2);
            else
                fprintf('%.2f\t%d\tfit_failed\tfit_failed\tfit_failed\tfit_failed\n', ...
                        ir_snr_points(snr_idx), n);
            end

            histogram(vals, edges_static, ...
                      'DisplayStyle', 'stairs', ...
                      'Normalization', 'probability', ...
                      'LineWidth', 1.6, ...
                      'EdgeColor', cmap_tx(i, :), ...
                      'DisplayName', sprintf('Tx n=%d (\\beta=%.2f)', n, beta_static_by_n(n+1)));
            has_data = true;
        end

        grid(ax, 'on');
        xlabel(ax, 'Effective SINR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(ax, 'Probability', 'FontSize', 12, 'FontWeight', 'bold');
        title(ax, sprintf('Input SNR = %.2f dB', ir_snr_points(snr_idx)), ...
              'FontSize', 12, 'FontWeight', 'bold');
        if has_data
            lgd = legend(ax, 'Location', 'best');
            lgd.FontSize = 10;
            lgd.FontWeight = 'bold';
        else
            text(ax, 0.5, 0.5, 'No samples', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
    end
    title(tl_static, 'Static \beta: Effective SINR Histogram by Input SNR (Tx n=0..3)', ...
          'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, 'effective_sinr_distribution_static_hist_tx0to3.png');
else
    warning('No valid static effective SINR samples found for Tx n=0..3 by-SNR histogram plot.');
end

vals_for_edges_dynamic = [];
for snr_idx = 1:n_snr_points
    for n = n_plot
        vals_n = snrEff_dynamic_by_n{n+1};
        snr_idx_n = sample_snr_idx_dynamic_by_n{n+1};
        if isempty(vals_n) || isempty(snr_idx_n)
            continue;
        end
        valid_mask = isfinite(vals_n) & (snr_idx_n == snr_idx);
        vals_for_edges_dynamic = [vals_for_edges_dynamic; vals_n(valid_mask)]; %#ok<AGROW>
    end
end

if ~isempty(vals_for_edges_dynamic)
    minEdgeDyn = floor(min(vals_for_edges_dynamic));
    maxEdgeDyn = ceil(max(vals_for_edges_dynamic));
    if minEdgeDyn == maxEdgeDyn
        minEdgeDyn = minEdgeDyn - 0.5;
        maxEdgeDyn = maxEdgeDyn + 0.5;
    end
    edges_dynamic = linspace(minEdgeDyn, maxEdgeDyn, 40);

    dyn_cols = min(3, n_snr_points);
    dyn_rows = ceil(n_snr_points / dyn_cols);
    figure('Name', 'Effective SINR Distribution: Dynamic Beta by SNR (Tx n=0..3)', ...
           'Position', [140, 140, 1400, 720]);
    tl = tiledlayout(dyn_rows, dyn_cols, 'TileSpacing', 'compact', 'Padding', 'compact');
    fprintf('\n--- Skew Generalized Normal Fits: Dynamic Beta ---\n');
    fprintf('SNR(dB)\tTx\tMean\t\tVariance\tLambda1\t\tLambda2\n');
    fprintf('---------------------------------------------------------------------\n');
    for snr_idx = 1:n_snr_points
        ax = nexttile(tl, snr_idx);
        hold(ax, 'on');
        has_data = false;
        for i = 1:numel(n_plot)
            n = n_plot(i);
            vals_n = snrEff_dynamic_by_n{n+1};
            snr_idx_n = sample_snr_idx_dynamic_by_n{n+1};
            if isempty(vals_n) || isempty(snr_idx_n)
                continue;
            end

            valid_mask = isfinite(vals_n) & (snr_idx_n == snr_idx);
            vals = vals_n(valid_mask);
            if isempty(vals)
                continue;
            end
            fit_dynamic = fitSkewGeneralizedNormal(vals);
            fit_dynamic_by_snr_tx{snr_idx, i} = fit_dynamic;
            if fit_dynamic.ok
                fprintf('%.2f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\n', ...
                        ir_snr_points(snr_idx), n, fit_dynamic.mean, fit_dynamic.variance, ...
                        fit_dynamic.lambda1, fit_dynamic.lambda2);
            else
                fprintf('%.2f\t%d\tfit_failed\tfit_failed\tfit_failed\tfit_failed\n', ...
                        ir_snr_points(snr_idx), n);
            end

            histogram(vals, edges_dynamic, ...
                      'DisplayStyle', 'stairs', ...
                      'Normalization', 'probability', ...
                      'LineWidth', 1.6, ...
                      'EdgeColor', cmap_tx(i, :), ...
                      'DisplayName', sprintf('Tx n=%d (\\beta=%.2f)', n, beta_dynamic_by_n(n+1)));
            has_data = true;
        end

        grid(ax, 'on');
        xlabel(ax, 'Effective SINR (dB)', 'FontSize', 12, 'FontWeight', 'bold');
        ylabel(ax, 'Probability', 'FontSize', 12, 'FontWeight', 'bold');
        title(ax, sprintf('Input SNR = %.2f dB', ir_snr_points(snr_idx)), ...
              'FontSize', 12, 'FontWeight', 'bold');
        if has_data
            lgd = legend(ax, 'Location', 'best');
            lgd.FontSize = 10;
            lgd.FontWeight = 'bold';
        else
            text(ax, 0.5, 0.5, 'No samples', ...
                 'Units', 'normalized', 'HorizontalAlignment', 'center');
        end
    end
    title(tl, 'Dynamic \beta: Effective SINR Histogram by Input SNR (Tx n=0..3)', ...
          'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, 'effective_sinr_distribution_dynamic_hist_tx0to3.png');
else
    warning('No valid dynamic effective SINR samples found for Tx n=0..3 by-SNR histogram plot.');
end

if ~isempty(vals_for_edges_dynamic)
    binw_dynamic = mean(diff(edges_dynamic));
    figure('Name', 'Skew Generalized Normal Fit Overlay (Dynamic Beta)', ...
           'Position', [170, 80, 1100, 1800]);
    tl_overlay = tiledlayout(5, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for snr_idx = 1:n_snr_points
        ax_d = nexttile(tl_overlay, snr_idx);
        hold(ax_d, 'on');
        has_data_d = false;
        for i = 1:numel(n_plot)
            n = n_plot(i);
            vals_n = snrEff_dynamic_by_n{n+1};
            snr_idx_n = sample_snr_idx_dynamic_by_n{n+1};
            if isempty(vals_n) || isempty(snr_idx_n)
                continue;
            end
            valid_mask = isfinite(vals_n) & (snr_idx_n == snr_idx);
            vals = vals_n(valid_mask);
            if isempty(vals)
                continue;
            end

            histogram(vals, edges_dynamic, ...
                      'DisplayStyle', 'stairs', ...
                      'Normalization', 'probability', ...
                      'LineStyle', '--', ...
                      'LineWidth', 1.2, ...
                      'EdgeColor', cmap_tx(i, :), ...
                      'HandleVisibility', 'off');

            fit_d = fit_dynamic_by_snr_tx{snr_idx, i};
            if ~isempty(fit_d) && isfield(fit_d, 'ok') && fit_d.ok
                xg = linspace(edges_dynamic(1), edges_dynamic(end), 500);
                yg = skewGenNormPdf(xg, fit_d.mu, fit_d.lambda1, fit_d.lambda2, fit_d.p) * binw_dynamic;
                legend_label = sprintf('Tx n=%d: mean=%.2f, var=%.2f, \\lambda_1=%.2f, \\lambda_2=%.2f', ...
                                       n, fit_d.mean, fit_d.variance, fit_d.lambda1, fit_d.lambda2);
                plot(ax_d, xg, yg, '-', 'Color', cmap_tx(i, :), ...
                     'LineWidth', 1.8, 'DisplayName', legend_label);
            else
                plot(ax_d, NaN, NaN, '-', 'Color', cmap_tx(i, :), ...
                     'LineWidth', 1.8, 'DisplayName', sprintf('Tx n=%d: fit failed', n));
            end
            has_data_d = true;
        end
        grid(ax_d, 'on');
        ax_d.FontSize = 12;
        ax_d.FontWeight = 'bold';
        ax_d.LineWidth = 1.2;
        xlabel(ax_d, 'Effective SINR (dB)', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel(ax_d, 'Probability', 'FontSize', 14, 'FontWeight', 'bold');
        title(ax_d, sprintf('Dynamic | SNR = %.2f dB', ir_snr_points(snr_idx)), ...
              'FontSize', 15, 'FontWeight', 'bold');
        if has_data_d
            lgd = legend(ax_d, 'Location', 'best');
            lgd.FontSize = 11;
            lgd.FontWeight = 'bold';
        else
            text(ax_d, 0.5, 0.5, 'No samples', 'Units', 'normalized', ...
                 'HorizontalAlignment', 'center');
        end
    end

    title(tl_overlay, 'Dynamic \beta: Histogram (dashed) and Skew Generalized Normal Fit (solid)', ...
          'FontSize', 18, 'FontWeight', 'bold');
    saveas(gcf, 'effective_sinr_distribution_fit_overlay.png');
else
    warning('Skipping fit-overlay figure because dynamic histogram data is incomplete.');
end

%% 5. SUMMARY TABLE
fprintf('\n================ FINAL SUMMARY ================\n');
fprintf('n\t| AWGN\t| beta_static\t| beta_dynamic\t| MSE_paper\t| MSE_dynamic\n');
fprintf('--------------------------------------------------------------------------\n');
for r = 1:numel(results)
    fprintf('%d\t| %s\t| %.4f\t\t| %.4f\t\t| %.6f\t| %.6f\n', ...
        results(r).n, results(r).awgn_ref, results(r).beta_static, ...
        results(r).beta_dynamic, results(r).mse_paper, results(r).mse_dynamic);
end
fprintf('==========================================================================\n');


% =========================================================================
% HELPER FUNCTIONS
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

function [snrVec, blerVec] = parseSummarySNRBler(summaryFile)
    fid = fopen(summaryFile, 'r');
    if fid < 0
        error('Cannot open file: %s', summaryFile);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

    lines = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
    lines = lines{1};

    snrVec = [];
    blerVec = [];
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
            blerVec(end+1,1) = vals(3); %#ok<AGROW>
        end
    end

    if isempty(snrVec)
        error('No SNR/BLER values found in %s', summaryFile);
    end

    [snrVec, ord] = sort(snrVec);
    blerVec = blerVec(ord);
    blerVec = min(max(blerVec, 0), 1);
end

function [retransInfo_all, row_snr_idx] = loadRetransInfoMatrix(folder)
    matFiles = dir(fullfile(folder, 'retrans_snr_*.mat'));
    if isempty(matFiles)
        error('No retrans_snr_*.mat files in %s', folder);
    end

    idxNum = inf(numel(matFiles), 1);
    for i = 1:numel(matFiles)
        tok = regexp(matFiles(i).name, 'retrans_snr_(\d+)\.mat', 'tokens', 'once');
        if ~isempty(tok)
            idxNum(i) = str2double(tok{1});
        end
    end
    [~, ord] = sort(idxNum);
    matFiles = matFiles(ord);
    idxNum = idxNum(ord);

    retransInfo_all = [];
    row_snr_idx = [];
    for i = 1:numel(matFiles)
        path_i = fullfile(folder, matFiles(i).name);
        d = load(path_i);
        if ~isfield(d, 'retransInfo')
            warning('Skipping %s (missing retransInfo)', path_i);
            continue;
        end
        nRows_i = size(d.retransInfo, 1);
        if isempty(retransInfo_all)
            retransInfo_all = d.retransInfo;
        else
            retransInfo_all = [retransInfo_all; d.retransInfo];
        end
        row_snr_idx = [row_snr_idx; repmat(idxNum(i), nRows_i, 1)]; %#ok<AGROW>
    end

    if isempty(retransInfo_all)
        error('No usable retransInfo found in %s', folder);
    end
end

function packetHistByRow = buildPacketHistories(retransInfo_all)
    nRows = size(retransInfo_all, 1);
    packetHistByRow = cell(nRows, 1);

    for i = 1:nRows
        row_data = retransInfo_all(i, :);
        packet_ids = [row_data.packetID];
        [~, ~, grp] = unique(packet_ids, 'stable');
        nPkt = max(grp);

        hist_i = repmat(struct('retransNum', [], 'idx', []), nPkt, 1);
        for g = 1:nPkt
            idx = find(grp == g);
            [~, ord] = sort([row_data(idx).retransNum]);
            idx = idx(ord);
            hist_i(g).idx = idx;
            hist_i(g).retransNum = [row_data(idx).retransNum];
        end
        packetHistByRow{i} = hist_i;
    end
end

function mse = awgnPerSnrFittingMse(sinrStore, perStore, lut, beta, useLogDomain, blerFloor)
    if beta <= 0 || beta > 1e6
        mse = 1e10;
        return;
    end

    [binsnr, binper, ~, ok] = fitPoints(sinrStore, perStore, beta);
    if ~ok
        mse = 1e10;
        return;
    end

    perIdeal = interp1(lut(:,1), lut(:,2), binsnr, 'linear', 'extrap');
    perIdeal = min(max(perIdeal(:), 0), 1);
    fitIdx = isfinite(binper) & (binper > 0);

    if sum(fitIdx) < 2
        mse = 1e10;
        return;
    end

    mse = computeBlerMismatch(binper(fitIdx), perIdeal(fitIdx), useLogDomain, blerFloor);
    if ~isfinite(mse)
        mse = 1e10;
    end
end

function snrEff = effectiveSinrVec(sinrStore, beta)
    if isvector(sinrStore), sinrStore = sinrStore(:)'; end
    linSinr = 10.^(sinrStore / 10);
    expTerm = exp(-linSinr / beta);
    expTerm(~isfinite(expTerm)) = NaN;
    meanExp = rowMeanOmitNan(expTerm);

    snrEffLin = -beta .* log(meanExp);
    snrEff = nan(size(snrEffLin));
    valid = isfinite(snrEffLin) & (snrEffLin > 0);
    snrEff(valid) = 10 * log10(snrEffLin(valid));
end

function valid = snrEffValidMask(sinrStore, perStore, beta)
    snrEffRaw = effectiveSinrVec(sinrStore, beta);
    valid = isfinite(snrEffRaw) & isfinite(perStore);
end

function fit = fitSkewGeneralizedNormal(x)
    x = x(:);
    x = x(isfinite(x));

    fit = struct('ok', false, ...
                 'n', numel(x), ...
                 'mu', NaN, ...
                 'p', NaN, ...
                 'lambda1', NaN, ...
                 'lambda2', NaN, ...
                 'mean', NaN, ...
                 'variance', NaN);

    if numel(x) < 6
        return;
    end

    mu0 = median(x);
    if ~isfinite(mu0)
        mu0 = mean(x);
    end

    std_all = std(x, 0);
    if ~isfinite(std_all) || std_all <= 0
        x_min = min(x);
        x_max = max(x);
        if isfinite(x_min) && isfinite(x_max) && (x_max > x_min)
            std_all = max(1e-3, (x_max - x_min) / 6);
        else
            std_all = 1.0;
        end
    end

    x_left = x(x <= mu0);
    x_right = x(x > mu0);
    lambda1_0 = std(x_left, 0);
    lambda2_0 = std(x_right, 0);
    if ~isfinite(lambda1_0) || lambda1_0 <= 0
        lambda1_0 = std_all;
    end
    if ~isfinite(lambda2_0) || lambda2_0 <= 0
        lambda2_0 = std_all;
    end

    theta0 = [mu0, log(lambda1_0), log(lambda2_0), log(2.0)];
    obj = @(th) skewGenNormNllFromTheta(th, x);
    opts = optimset('Display', 'off', ...
                    'MaxFunEvals', 800, ...
                    'MaxIter', 800, ...
                    'TolX', 1e-6, ...
                    'TolFun', 1e-6);

    [theta_hat, fval, exitflag] = fminsearch(obj, theta0, opts);
    if exitflag <= 0 || ~isfinite(fval)
        theta_hat = theta0;
    end

    [nll, params] = skewGenNormNllFromTheta(theta_hat, x);
    if ~isfinite(nll) || ~params.valid
        return;
    end

    [mean_fit, var_fit] = skewGenNormMoments(params.mu, params.lambda1, params.lambda2, params.p);
    if ~isfinite(mean_fit) || ~isfinite(var_fit)
        return;
    end

    fit.ok = true;
    fit.mu = params.mu;
    fit.p = params.p;
    fit.lambda1 = params.lambda1;
    fit.lambda2 = params.lambda2;
    fit.mean = mean_fit;
    fit.variance = max(var_fit, 0);
end

function [nll, params] = skewGenNormNllFromTheta(theta, x)
    mu = theta(1);
    lambda1 = exp(theta(2));
    lambda2 = exp(theta(3));
    p = exp(theta(4));

    params = struct('valid', false, 'mu', mu, 'lambda1', lambda1, 'lambda2', lambda2, 'p', p);
    nll_bad = 1e12;

    if ~isfinite(mu) || ~isfinite(lambda1) || ~isfinite(lambda2) || ~isfinite(p)
        nll = nll_bad;
        return;
    end
    if lambda1 <= 1e-9 || lambda2 <= 1e-9 || p <= 0.2 || p >= 20
        nll = nll_bad + 1e6 * (abs(log(max(p, eps))) + 1);
        return;
    end

    lambda_x = lambda2 * ones(size(x));
    lambda_x(x <= mu) = lambda1;
    z = abs(x - mu) ./ lambda_x;
    if any(~isfinite(z))
        nll = nll_bad;
        return;
    end

    log_c = log(p) - log(lambda1 + lambda2) - gammaln(1 / p);
    log_pdf = log_c - z.^p;
    if any(~isfinite(log_pdf))
        nll = nll_bad;
        return;
    end

    nll = -sum(log_pdf);
    if ~isfinite(nll)
        nll = nll_bad;
        return;
    end

    params.valid = true;
end

function [mean_fit, var_fit] = skewGenNormMoments(mu, lambda1, lambda2, p)
    r21 = exp(gammaln(2 / p) - gammaln(1 / p));
    r31 = exp(gammaln(3 / p) - gammaln(1 / p));

    delta = (lambda2 - lambda1) * r21;
    mean_fit = mu + delta;
    second_about_mu = ((lambda1^3 + lambda2^3) / (lambda1 + lambda2)) * r31;
    var_fit = second_about_mu - delta^2;
end

function y = skewGenNormPdf(x, mu, lambda1, lambda2, p)
    y = zeros(size(x));
    if ~isfinite(mu) || ~isfinite(lambda1) || ~isfinite(lambda2) || ~isfinite(p)
        y(:) = NaN;
        return;
    end
    if lambda1 <= 0 || lambda2 <= 0 || p <= 0
        y(:) = NaN;
        return;
    end

    lam = lambda2 * ones(size(x));
    lam(x <= mu) = lambda1;
    z = abs(x - mu) ./ lam;
    c = p / ((lambda1 + lambda2) * gamma(1 / p));
    y = c * exp(-(z.^p));
end

function M = cellsToMatrix(cellsIn)
    n = numel(cellsIn);
    lens = cellfun(@numel, cellsIn);
    maxLen = max(lens);
    minLen = min(lens);

    if minLen == maxLen
        M = vertcat(cellsIn{:});
        return;
    end

    M = NaN(n, maxLen);
    for i = 1:n
        li = lens(i);
        if li > 0
            M(i, 1:li) = cellsIn{i};
        end
    end
end

function [sinrOut, perOut] = downsamplePackets(sinrIn, perIn, maxN)
    n = size(sinrIn, 1);
    if n <= maxN
        sinrOut = sinrIn;
        perOut = perIn;
        return;
    end

    idx = randperm(n, maxN);
    sinrOut = sinrIn(idx, :);
    perOut = perIn(idx, :);
end

function y = rowMeanOmitNan(x)
    mask = ~isnan(x);
    cnt = sum(mask, 2);
    x(~mask) = 0;
    y = sum(x, 2) ./ cnt;
    y(cnt == 0) = NaN;
end

function [binsnr, binper, snrEff, mse] = awgnPerSnrFitting(sinrStore, perStore, lut, beta, useLogDomain, blerFloor)
    [binsnr, binper, snrEff, ok] = fitPoints(sinrStore, perStore, beta);
    if ~ok
        mse = 1e10;
        return;
    end

    perIdeal = interp1(lut(:,1), lut(:,2), binsnr, 'linear', 'extrap');
    perIdeal = min(max(perIdeal(:), 0), 1);
    fitIdx = isfinite(binper) & (binper > 0);
    if sum(fitIdx) < 2
        mse = 1e10;
    else
        mse = computeBlerMismatch(binper(fitIdx), perIdeal(fitIdx), useLogDomain, blerFloor);
    end
end

function mse = computeBlerMismatch(blerMeas, blerRef, useLogDomain, blerFloor)
    if nargin < 4 || isempty(blerFloor)
        blerFloor = 1e-5;
    end
    blerMeas = min(max(blerMeas(:), blerFloor), 1);
    blerRef = min(max(blerRef(:), blerFloor), 1);

    if useLogDomain
        d = log10(blerMeas) - log10(blerRef);
    else
        d = blerMeas - blerRef;
    end
    mse = mean(d.^2);
end

function [binsnr, binper, snrEff, ok] = fitPoints(sinrStore, perStore, beta)
    snrEff = effectiveSinrVec(sinrStore, beta);
    valid = isfinite(snrEff) & isfinite(perStore);
    snrEff = snrEff(valid);
    perStore = perStore(valid);

    binsnr = [];
    binper = [];
    ok = false;

    if numel(snrEff) < 2
        return;
    end

    edges = floor(min(snrEff)):0.25:ceil(max(snrEff));
    if numel(edges) < 2
        return;
    end

    binsnr = edges(1:end-1) + (edges(2)-edges(1))/2;
    BINS = discretize(snrEff, edges);
    binper = nan(numel(edges)-1, 1);
    for i = 1:numel(edges)-1
        m = (BINS == i);
        c = sum(m);
        if c > 0
            binper(i) = mean(perStore(m));
        end
    end
    ok = true;
end
