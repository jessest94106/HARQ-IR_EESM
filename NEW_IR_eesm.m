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
manual_ir_folder = '';

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
retransInfo_ir_all = loadRetransInfoMatrix(ir_folder);
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
        end
    end

    if n_keep == 0
        fprintf('No IR samples for n = %d. Skipping.\n', target_n);
        continue;
    end
    sinr_cells = sinr_cells(1:n_keep);
    perStore_n = perStore_n(1:n_keep);
    sinrStore_n = cellsToMatrix(sinr_cells);

    [sinrFit_n, perFit_n] = downsamplePackets(sinrStore_n, perStore_n, max_fit_packets);
    fprintf('Samples: full=%d, fit=%d\n', size(sinrStore_n,1), size(sinrFit_n,1));

    subplot(rows, cols, target_n + 1);

    if target_n == 0
        mse_func = @(beta) awgnPerSnrFittingMse(sinrFit_n, perFit_n, lut, beta, use_log_bler_objective, bler_floor);
        beta_static = fminbnd(mse_func, beta_lb, beta_ub, options);
        [binsnr, binper, ~, mse_val] = awgnPerSnrFitting(sinrStore_n, perStore_n, lut, beta_static, use_log_bler_objective, bler_floor);

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
        [binsnr_p, binper_p, ~, mse_p] = awgnPerSnrFitting(sinrStore_n, perStore_n, lut, beta_static, use_log_bler_objective, bler_floor);

        mse_func_dyn = @(beta) awgnPerSnrFittingMse(sinrFit_n, perFit_n, lut, beta, use_log_bler_objective, bler_floor);
        beta_dyn = fminbnd(mse_func_dyn, beta_lb, beta_ub, options);
        [binsnr_d, binper_d, ~, mse_d] = awgnPerSnrFitting(sinrStore_n, perStore_n, lut, beta_dyn, use_log_bler_objective, bler_floor);

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

%% 4. SUMMARY TABLE
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

function retransInfo_all = loadRetransInfoMatrix(folder)
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

    retransInfo_all = [];
    for i = 1:numel(matFiles)
        path_i = fullfile(folder, matFiles(i).name);
        d = load(path_i);
        if ~isfield(d, 'retransInfo')
            warning('Skipping %s (missing retransInfo)', path_i);
            continue;
        end
        if isempty(retransInfo_all)
            retransInfo_all = d.retransInfo;
        else
            retransInfo_all = [retransInfo_all; d.retransInfo];
        end
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
