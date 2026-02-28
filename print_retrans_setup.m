function print_retrans_setup(pattern)
%PRINT_RETRANS_SETUP Print setup and retransmission stats from MAT files.
%   PRINT_RETRANS_SETUP() looks for ./arch_results/*retrans_snr_1.mat first,
%   then falls back to common project layouts under ./result/ and ./results/.
%
%   PRINT_RETRANS_SETUP(PATTERN) uses the provided dir() pattern.

if nargin < 1 || isempty(pattern)
    pattern = './arch_results/*retrans_snr_1.mat';
end

fileList = iResolveMatFiles(pattern);
if isempty(fileList)
    error('No files matched pattern "%s" (or fallback patterns).', pattern);
end

fprintf('Found %d file(s).\n', numel(fileList));
for idx = 1:numel(fileList)
    filePath = fileList{idx};
    runFolder = fileparts(filePath);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s\n', idx, numel(fileList), filePath);
    fprintf('============================================================\n');

    data = load(filePath);
    retransInfo = [];
    if isfield(data, 'retransInfo') && ~isempty(data.retransInfo)
        retransInfo = data.retransInfo;
    else
        fprintf('No retransInfo payload in this file.\n');
    end

    setup = iGetSetup(filePath);
    iPrintSetup(setup);
    if ~isempty(retransInfo)
        iPrintRetransStats(retransInfo);
    end

    settingPath = fullfile(runFolder, 'setting.txt');
    iWriteSettingFile(settingPath, setup, retransInfo);
    fprintf('Saved: %s\n', settingPath);
end
end

function fileList = iResolveMatFiles(primaryPattern)
patterns = { ...
    primaryPattern, ...
    './arch_results/*/retrans_snr_1.mat', ...
    './arch_results/**/retrans_snr_1.mat' ...
    };

fileList = {};
for i = 1:numel(patterns)
    d = dir(patterns{i});
    d = d(~[d.isdir]);
    if ~isempty(d)
        fileList = cell(numel(d), 1);
        for k = 1:numel(d)
            fileList{k} = fullfile(d(k).folder, d(k).name);
        end
        return;
    end
end
end

function setup = iGetSetup(matFilePath)
setup = struct( ...
    'DelayProfile', '', ...
    'DelaySpread', NaN, ...
    'MaximumDopplerShift', NaN, ...
    'NTxAnts', NaN, ...
    'NRxAnts', NaN, ...
    'NFrames', NaN, ...
    'TotalSlots', NaN, ...
    'Modulation', '', ...
    'HARQ', '' ...
    );

summaryPath = fullfile(fileparts(matFilePath), 'simulation_summary.txt');
if exist(summaryPath, 'file')
    setup = iParseSummary(summaryPath, setup);
end

% Fill any missing setup items from folder name if needed.
parentFolder = fileparts(matFilePath);
[~, runName] = fileparts(parentFolder);
setup = iBackfillFromFolderName(setup, runName);
end

function setup = iParseSummary(summaryPath, setup)
fid = fopen(summaryPath, 'r');
if fid < 0
    return;
end

cleanup = onCleanup(@() fclose(fid));
while true
    line = fgetl(fid);
    if ~ischar(line)
        break;
    end

    token = regexp(line, '^Delay Profile:\s*(.+)$', 'tokens', 'once');
    if ~isempty(token)
        setup.DelayProfile = strtrim(token{1});
        continue;
    end

    token = regexp(line, '^Delay Spread:\s*([-\deE.+]+)', 'tokens', 'once');
    if ~isempty(token)
        setup.DelaySpread = str2double(token{1});
        continue;
    end

    token = regexp(line, '^Maximum Doppler Shift:\s*([-\deE.+]+)', 'tokens', 'once');
    if ~isempty(token)
        setup.MaximumDopplerShift = str2double(token{1});
        continue;
    end

    token = regexp(line, '^NTxAnts:\s*(\d+),\s*NRxAnts:\s*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        setup.NTxAnts = str2double(token{1});
        setup.NRxAnts = str2double(token{2});
        continue;
    end

    token = regexp(line, '^NFrames:\s*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        setup.NFrames = str2double(token{1});
        continue;
    end

    token = regexp(line, '^Total Slots:\s*(\d+)', 'tokens', 'once');
    if ~isempty(token)
        setup.TotalSlots = str2double(token{1});
        continue;
    end

    token = regexp(line, '^Modulation:\s*(.+)$', 'tokens', 'once');
    if ~isempty(token)
        setup.Modulation = strtrim(token{1});
        continue;
    end

    token = regexp(line, '^HARQ:\s*(.+)$', 'tokens', 'once');
    if ~isempty(token)
        setup.HARQ = strtrim(token{1});
    end
end
end

function setup = iBackfillFromFolderName(setup, runName)
token = regexp(runName, '^.*_(?<harq>[^_]+)_(?<mod>[^_]+(?:_[^_]+)*)_(?<tx>\d+)x(?<rx>\d+)_(?<ds>[^_]+)$', 'names', 'once');
if isempty(token)
    return;
end

if isempty(setup.HARQ)
    setup.HARQ = token.harq;
end

if isempty(setup.Modulation)
    setup.Modulation = strrep(token.mod, '_', '/');
end

if isnan(setup.NTxAnts)
    setup.NTxAnts = str2double(token.tx);
end

if isnan(setup.NRxAnts)
    setup.NRxAnts = str2double(token.rx);
end

if isnan(setup.DelaySpread)
    setup.DelaySpread = str2double(token.ds);
end
end

function iPrintSetup(setup)
fprintf('Setup\n');
fprintf('  Antenna setup: %s\n', iAntennaString(setup.NTxAnts, setup.NRxAnts));
fprintf('  HARQ mode: %s\n', iValOrUnknown(setup.HARQ));
fprintf('  Modulation: %s\n', iValOrUnknown(setup.Modulation));
fprintf('  Delay profile: %s\n', iValOrUnknown(setup.DelayProfile));

if isnan(setup.DelaySpread)
    fprintf('  Delay spread: unknown\n');
else
    fprintf('  Delay spread: %.3e s\n', setup.DelaySpread);
end

if isnan(setup.MaximumDopplerShift)
    fprintf('  Max Doppler shift: unknown\n');
else
    fprintf('  Max Doppler shift: %.2f Hz\n', setup.MaximumDopplerShift);
end

if isnan(setup.NFrames)
    fprintf('  NFrames: unknown\n');
else
    fprintf('  NFrames: %d\n', setup.NFrames);
end

if isnan(setup.TotalSlots)
    fprintf('  Total slots: unknown\n');
else
    fprintf('  Total slots: %d\n', setup.TotalSlots);
end
end

function iPrintRetransStats(retransInfo)
packetID = [retransInfo.packetID];
retransNum = [retransInfo.retransNum];
errorState = [retransInfo.errorState];

fprintf('Retransmission Stats\n');
fprintf('  Records: %d\n', numel(retransInfo));
fprintf('  Unique packets: %d\n', numel(unique(packetID)));
fprintf('  Retransmission number range: %d to %d\n', min(retransNum), max(retransNum));

uRetrans = unique(retransNum);
for i = 1:numel(uRetrans)
    r = uRetrans(i);
    mask = retransNum == r;
    errRate = mean(errorState(mask));
    fprintf('  Retrans #%d: count=%d, error-rate=%.4f\n', r, sum(mask), errRate);
end

% Final packet outcome: take the highest retransmission index per packet.
uPkt = unique(packetID);
finalErr = zeros(numel(uPkt), 1);
for i = 1:numel(uPkt)
    idx = find(packetID == uPkt(i));
    [~, localMaxPos] = max(retransNum(idx));
    finalErr(i) = errorState(idx(localMaxPos));
end
fprintf('  Final packet BLER: %.4f\n', mean(finalErr));

if isfield(retransInfo, 'sinr')
    sinr = [retransInfo.sinr];
    fprintf('  SINR (dB): mean=%.3f, min=%.3f, max=%.3f\n', mean(sinr), min(sinr), max(sinr));
end
end

function out = iAntennaString(nTx, nRx)
if isnan(nTx) || isnan(nRx)
    out = 'unknown';
else
    out = sprintf('%dx%d (NTx x NRx)', nTx, nRx);
end
end

function out = iValOrUnknown(val)
if isempty(val)
    out = 'unknown';
else
    out = val;
end
end

function iWriteSettingFile(settingPath, setup, retransInfo)
fid = fopen(settingPath, 'w');
if fid < 0
    warning('Could not create %s', settingPath);
    return;
end

cleanup = onCleanup(@() fclose(fid));
fprintf(fid, 'Setup\n');
fprintf(fid, 'Antenna setup: %s\n', iAntennaString(setup.NTxAnts, setup.NRxAnts));
fprintf(fid, 'HARQ mode: %s\n', iValOrUnknown(setup.HARQ));
fprintf(fid, 'Modulation: %s\n', iValOrUnknown(setup.Modulation));
fprintf(fid, 'Delay profile: %s\n', iValOrUnknown(setup.DelayProfile));

if isnan(setup.DelaySpread)
    fprintf(fid, 'Delay spread: unknown\n');
else
    fprintf(fid, 'Delay spread: %.3e s\n', setup.DelaySpread);
end

if isnan(setup.MaximumDopplerShift)
    fprintf(fid, 'Max Doppler shift: unknown\n');
else
    fprintf(fid, 'Max Doppler shift: %.2f Hz\n', setup.MaximumDopplerShift);
end

if isnan(setup.NFrames)
    fprintf(fid, 'NFrames: unknown\n');
else
    fprintf(fid, 'NFrames: %d\n', setup.NFrames);
end

if isnan(setup.TotalSlots)
    fprintf(fid, 'Total slots: unknown\n');
else
    fprintf(fid, 'Total slots: %d\n', setup.TotalSlots);
end

if isempty(retransInfo)
    fprintf(fid, '\nRetransmission stats: unavailable\n');
    return;
end

packetID = [retransInfo.packetID];
retransNum = [retransInfo.retransNum];
errorState = [retransInfo.errorState];

fprintf(fid, '\nRetransmission Stats\n');
fprintf(fid, 'Records: %d\n', numel(retransInfo));
fprintf(fid, 'Unique packets: %d\n', numel(unique(packetID)));
fprintf(fid, 'Retransmission number range: %d to %d\n', min(retransNum), max(retransNum));

uRetrans = unique(retransNum);
for i = 1:numel(uRetrans)
    r = uRetrans(i);
    mask = retransNum == r;
    errRate = mean(errorState(mask));
    fprintf(fid, 'Retrans #%d: count=%d, error-rate=%.4f\n', r, sum(mask), errRate);
end

uPkt = unique(packetID);
finalErr = zeros(numel(uPkt), 1);
for i = 1:numel(uPkt)
    idx = find(packetID == uPkt(i));
    [~, localMaxPos] = max(retransNum(idx));
    finalErr(i) = errorState(idx(localMaxPos));
end
fprintf(fid, 'Final packet BLER: %.4f\n', mean(finalErr));

if isfield(retransInfo, 'sinr')
    sinr = [retransInfo.sinr];
    fprintf(fid, 'SINR (dB): mean=%.3f, min=%.3f, max=%.3f\n', mean(sinr), min(sinr), max(sinr));
end
end
