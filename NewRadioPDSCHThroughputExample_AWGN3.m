%% NR PDSCH Throughput 

simParameters = struct();       % Clear simParameters variable to contain all key simulation parameters 
simParameters.NFrames = 1000;      % Number of 10 ms frames
% simParameters.SNRIn = [linspace(-10.5,6.5, 16) 6.8 7.8 8.5 8.9 9.3]; % SNR range (dB)
simParameters.SNRIn = [linspace(6.5,7.62, 5)]; % SNR range (dB)

%% Channel Estimator Configuration
% The logical variable |PerfectChannelEstimator| controls channel
% estimation and synchronization behavior. When set to |true|, perfect
% channel estimation and synchronization is used. Otherwise, practical
% channel estimation and synchronization is used, based on the values of
% the received PDSCH DM-RS.

simParameters.PerfectChannelEstimator = true;

%% Simulation Diagnostics
% The variable |DisplaySimulationInformation| controls the display of
% simulation information such as the HARQ process ID used for each
% subframe. In case of CRC error, the value of the index to the RV sequence
% is also displayed.

simParameters.DisplaySimulationInformation = true;

simParameters.DisplayDiagnostics = false;

%% Carrier and PDSCH Configuration
% Set the key parameters of the simulation. These include:
% 
% * The bandwidth in resource blocks (12 subcarriers per resource block).
% * Subcarrier spacing: 15, 30, 60, 120 (kHz)
% * Cyclic prefix length: normal or extended
% * Cell ID
% * Number of transmit and receive antennas
% 
% A substructure containing the DL-SCH and PDSCH parameters is also
% specified. This includes:
% 
% * Target code rate
% * Allocated resource blocks (PRBSet)
% * Modulation scheme: 'QPSK', '16QAM', '64QAM', '256QAM'
% * Number of layers
% * PDSCH mapping type
% * DM-RS configuration parameters
% * PT-RS configuration parameters
% 
% Other simulation wide parameters are:
% 
% * Propagation channel model delay profile (TDL or CDL)

% Set waveform type and PDSCH numerology (SCS and CP type)
simParameters.Carrier = nrCarrierConfig;         % Carrier resource grid configuration
simParameters.Carrier.NSizeGrid = 51;            % Bandwidth in number of resource blocks (51 RBs at 30 kHz SCS for 20 MHz BW)
simParameters.Carrier.SubcarrierSpacing = 30;    % 15, 30, 60, 120 (kHz)
simParameters.Carrier.CyclicPrefix = 'Normal';   % 'Normal' or 'Extended' (Extended CP is relevant for 60 kHz SCS only)
simParameters.Carrier.NCellID = 1;               % Cell identity

% PDSCH/DL-SCH parameters
simParameters.PDSCH = nrPDSCHConfig;      % This PDSCH definition is the basis for all PDSCH transmissions in the BLER simulation
simParameters.PDSCHExtension = struct();  % This structure is to hold additional simulation parameters for the DL-SCH and PDSCH

% Define PDSCH time-frequency resource allocation per slot to be full grid (single full grid BWP)
simParameters.PDSCH.PRBSet = 0:simParameters.Carrier.NSizeGrid-1;                 % PDSCH PRB allocation
simParameters.PDSCH.SymbolAllocation = [0,simParameters.Carrier.SymbolsPerSlot];  % Starting symbol and number of symbols of each PDSCH allocation
simParameters.PDSCH.MappingType = 'A';     % PDSCH mapping type ('A'(slot-wise),'B'(non slot-wise))

% Scrambling identifiers
simParameters.PDSCH.NID = simParameters.Carrier.NCellID;
simParameters.PDSCH.RNTI = 1;

% PDSCH resource block mapping (TS 38.211 Section 7.3.1.6)
simParameters.PDSCH.VRBToPRBInterleaving = 0; % Disable interleaved resource mapping
simParameters.PDSCH.VRBBundleSize = 4;

% Define the number of transmission layers to be used
simParameters.PDSCH.NumLayers = 1;            % Number of PDSCH transmission layers

% Define codeword modulation and target coding rate
% The number of codewords is directly dependent on the number of layers so ensure that 
% layers are set first before getting the codeword number
if simParameters.PDSCH.NumCodewords > 1                             % Multicodeword transmission (when number of layers being > 4)
    simParameters.PDSCH.Modulation = {'16QAM','16QAM'};             % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PDSCHExtension.TargetCodeRate = [490 490]/1024;   % Code rate used to calculate transport block sizes
else
    simParameters.PDSCH.Modulation = '64QAM';                       % 'QPSK', '16QAM', '64QAM', '256QAM'
    simParameters.PDSCHExtension.TargetCodeRate = (772/1024)/3;         % Code rate used to calculate transport block sizes
end

% DM-RS and antenna port configuration (TS 38.211 Section 7.4.1.1)
simParameters.PDSCH.DMRS.DMRSPortSet = 0:simParameters.PDSCH.NumLayers-1; % DM-RS ports to use for the layers
simParameters.PDSCH.DMRS.DMRSTypeAPosition = 2;      % Mapping type A only. First DM-RS symbol position (2,3)
simParameters.PDSCH.DMRS.DMRSLength = 1;             % Number of front-loaded DM-RS symbols (1(single symbol),2(double symbol))
simParameters.PDSCH.DMRS.DMRSAdditionalPosition = 2; % Additional DM-RS symbol positions (max range 0...3)
simParameters.PDSCH.DMRS.DMRSConfigurationType = 2;  % DM-RS configuration type (1,2)
simParameters.PDSCH.DMRS.NumCDMGroupsWithoutData = 1;% Number of CDM groups without data
simParameters.PDSCH.DMRS.NIDNSCID = 1;               % Scrambling identity (0...65535)
simParameters.PDSCH.DMRS.NSCID = 0;                  % Scrambling initialization (0,1)

% PT-RS configuration (TS 38.211 Section 7.4.1.2)
simParameters.PDSCH.EnablePTRS = 0;                  % Enable or disable PT-RS (1 or 0)
simParameters.PDSCH.PTRS.TimeDensity = 1;            % PT-RS time density (L_PT-RS) (1, 2, 4)
simParameters.PDSCH.PTRS.FrequencyDensity = 2;       % PT-RS frequency density (K_PT-RS) (2 or 4)
simParameters.PDSCH.PTRS.REOffset = '00';            % PT-RS resource element offset ('00', '01', '10', '11')
simParameters.PDSCH.PTRS.PTRSPortSet = [];           % PT-RS antenna port, subset of DM-RS port set. Empty corresponds to lower DM-RS port number

% Reserved PRB patterns, if required (for CORESETs, forward compatibility etc)
simParameters.PDSCH.ReservedPRB{1}.SymbolSet = [];   % Reserved PDSCH symbols
simParameters.PDSCH.ReservedPRB{1}.PRBSet = [];      % Reserved PDSCH PRBs
simParameters.PDSCH.ReservedPRB{1}.Period = [];      % Periodicity of reserved resources

% Additional simulation and DL-SCH related parameters
%
% PDSCH PRB bundling (TS 38.214 Section 5.1.2.3)
simParameters.PDSCHExtension.PRGBundleSize = [];     % 2, 4, or [] to signify "wideband"
%
% HARQ process and rate matching/TBS parameters
simParameters.PDSCHExtension.XOverhead = 6*simParameters.PDSCH.EnablePTRS; % Set PDSCH rate matching overhead for TBS (Xoh) to 6 when PT-RS is enabled, otherwise 0
simParameters.PDSCHExtension.NHARQProcesses = 16;    % Number of parallel HARQ processes to use
simParameters.PDSCHExtension.EnableHARQ = "nan";      % Enable retransmissions for each process, using RV sequence [0,2,3,1]
if simParameters.PDSCHExtension.EnableHARQ=="IR"
    % In the final report of RAN WG1 meeting #91 (R1-1719301), it was
    % observed in R1-1717405 that if performance is the priority, [0 2 3 1]
    % should be used. If self-decodability is the priority, it should be
    % taken into account that the upper limit of the code rate at which
    % each RV is self-decodable is in the following order: 0>3>2>1
    rvSeq = [0 2 3 1];
    % rvSeq = [0 2];

elseif simParameters.PDSCHExtension.EnableHARQ=="CC"
    rvSeq = [0 0]; % Chase combining - all transmissions use RV=0
else
    % HARQ disabled - single transmission with RV=0, no retransmissions
    rvSeq = 0;
end
% LDPC decoder parameters
% Available algorithms: 'Belief propagation', 'Layered belief propagation', 'Normalized min-sum', 'Offset min-sum'
simParameters.PDSCHExtension.LDPCDecodingAlgorithm = 'Normalized min-sum';
simParameters.PDSCHExtension.MaximumLDPCIterationCount = 6;

% Define the overall transmission antenna geometry at end-points
% If using a CDL propagation channel then the integer number of antenna elements is
% turned into an antenna panel configured when the channel model object is created
simParameters.NTxAnts = 1;                        % Number of PDSCH transmission antennas (1,2,4,8,16,32,64,128,256,512,1024) >= NumLayers
if simParameters.PDSCH.NumCodewords > 1           % Multi-codeword transmission
    simParameters.NRxAnts = 8;                    % Number of UE receive antennas (even number >= NumLayers)
else
    simParameters.NRxAnts = 1;                    % Number of UE receive antennas (1 or even number >= NumLayers)
end

% Define data type ('single' or 'double') for resource grids and waveforms
simParameters.DataType = 'single';

% Define the general CDL/TDL propagation channel parameters
simParameters.DelayProfile = 'TDL-A';   % Use CDL-C model (Urban macrocell model)
simParameters.DelaySpread = 0;   % Channel delay spread (seconds)
simParameters.MaximumDopplerShift = 0; % Maximum Doppler shift (Hz)

% Cross-check the PDSCH layering against the channel geometry 
validateNumLayers(simParameters);

%% 
% The simulation relies on various pieces of information about the baseband 
% waveform, such as sample rate.

waveformInfo = nrOFDMInfo(simParameters.Carrier); % Get information about the baseband waveform after OFDM modulation step

%% Propagation Channel Model Construction
% Create the channel model object for the simulation. Both CDL and TDL channel 
% models are supported [ <#13 5> ].

% Constructed the CDL or TDL channel model object
if contains(simParameters.DelayProfile,'CDL','IgnoreCase',true)
    
    channel = nrCDLChannel; % CDL channel object
    
    % Turn the number of antennas into antenna panel array layouts. If
    % NTxAnts is not one of (1,2,4,8,16,32,64,128,256,512,1024), its value
    % is rounded up to the nearest value in the set. If NRxAnts is not 1 or
    % even, its value is rounded up to the nearest even number.
    channel = hArrayGeometry(channel,simParameters.NTxAnts,simParameters.NRxAnts);
    simParameters.NTxAnts = prod(channel.TransmitAntennaArray.Size);
    simParameters.NRxAnts = prod(channel.ReceiveAntennaArray.Size);
else
    channel = nrTDLChannel; % TDL channel object
    
    % Set the channel geometry
    channel.NumTransmitAntennas = simParameters.NTxAnts;
    channel.NumReceiveAntennas = simParameters.NRxAnts;
end

% Assign simulation channel parameters and waveform sample rate to the object
channel.DelayProfile = simParameters.DelayProfile;
channel.DelaySpread = simParameters.DelaySpread;
channel.MaximumDopplerShift = simParameters.MaximumDopplerShift;
channel.SampleRate = waveformInfo.SampleRate;

%% 
% Get the maximum number of delayed samples by a channel multipath
% component. This is calculated from the channel path with the largest
% delay and the implementation delay of the channel filter. This is
% required later to flush the channel filter to obtain the received signal.

chInfo = info(channel);
maxChDelay = chInfo.MaximumChannelDelay;

%% Processing Loop
% To determine the throughput at each SNR point, analyze the PDSCH data 
% per transmission instance using the following steps:
% 
% * _Update current HARQ process._ Check the transmission status for the
% given HARQ process to determine whether a retransmission is required. If
% that is not the case then generate new data.
% * _Resource grid generation._ Perform channel coding by calling the
% <docid:5g_ref#mw_sysobj_nrDLSCH nrDLSCH> System object. The object
% operates on the input transport block and keeps an internal copy of the
% transport block in case a retransmission is required. Modulate the coded
% bits on the PDSCH by using the <docid:5g_ref#mw_function_nrPDSCH nrPDSCH>
% function. Then apply precoding to the resulting signal.
% * _Waveform generation._ OFDM modulate the generated grid.
% * _Noisy channel modeling._ Pass the waveform through a CDL or TDL
% fading channel. Add AWGN. For an explanation of the SNR definition that
% this example uses, see
% <docid:5g_ug#mw_37cef3ca-2f4b-433d-8d68-117a881ca5fd SNR Definition used
% in Link Simulations>.
% * _Perform synchronization and OFDM demodulation._ For perfect
% synchronization, reconstruct the channel impulse response to synchronize
% the received waveform. For practical synchronization, correlate the
% received waveform with the PDSCH DM-RS. Then OFDM demodulate the
% synchronized signal.
% * _Perform channel estimation._ For perfect channel estimation,
% reconstruct the channel impulse response and perform OFDM demodulation.
% For practical channel estimation, use the PDSCH DM-RS.
% * _Perform equalization and CPE compensation._ MMSE equalize the
% estimated channel. Estimate the common phase error (CPE) by using the
% PT-RS symbols, then correct the error in each OFDM symbol within the
% range of reference PT-RS OFDM symbols.
% * _Precoding matrix calculation._ Generate the precoding matrix W for the
% next transmission by using singular value decomposition (SVD).
% * _Decode the PDSCH._ To obtain an estimate of the received codewords,
% demodulate and descramble the recovered PDSCH symbols for all transmit
% and receive antenna pairs, along with a noise estimate, by using the
% <docid:5g_ref#mw_function_nrPDSCHDecode nrPDSCHDecode> function.
% * _Decode the downlink shared channel (DL-SCH) and update HARQ process
% with the block CRC error._ Pass the vector of decoded soft bits to the
% <docid:5g_ref#mw_sysobj_nrDLSCHDecoder nrDLSCHDecoder> System object. The
% object decodes the codeword and returns the block CRC error used to
% determine the throughput of the system.

% Array to store the maximum throughput for all SNR points
maxThroughput = zeros(length(simParameters.SNRIn),1); 
% Array to store the simulation throughput for all SNR points
simThroughput = zeros(length(simParameters.SNRIn),1);


% Create DL-SCH encoder system object to perform transport channel encoding
encodeDLSCH = nrDLSCH;
encodeDLSCH.MultipleHARQProcesses = true;
encodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;

% Create DL-SCH decoder system object to perform transport channel decoding
% Use layered belief propagation for LDPC decoding, with half the number of
% iterations as compared to the default for belief propagation decoding
decodeDLSCH = nrDLSCHDecoder;
decodeDLSCH.MultipleHARQProcesses = true;
decodeDLSCH.TargetCodeRate = simParameters.PDSCHExtension.TargetCodeRate;
decodeDLSCH.LDPCDecodingAlgorithm = simParameters.PDSCHExtension.LDPCDecodingAlgorithm;
decodeDLSCH.MaximumLDPCIterationCount = simParameters.PDSCHExtension.MaximumLDPCIterationCount;

% Create a subfolder for saving plots
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');
harqMode = simParameters.PDSCHExtension.EnableHARQ;
modScheme = strjoin(string(simParameters.PDSCH.Modulation), '_');
subfolderName = sprintf('%s_%s_%s_%dx%d_%d', timestamp, harqMode, modScheme, simParameters.NTxAnts, simParameters.NRxAnts, simParameters.DelaySpread);
resultsDir = fullfile('AWGN_results', subfolderName);
if ~exist('AWGN_results', 'dir')
    mkdir('AWGN_results');  
end
if ~exist(resultsDir, 'dir')
    mkdir(resultsDir);
end

% Initialize cell array to collect retransmission data from each SNR point
retransInfoAll = cell(numel(simParameters.SNRIn), 1);

parfor snrIdx = 1:numel(simParameters.SNRIn)      % comment out for parallel computing
% parfor snrIdx = 1:numel(simParameters.SNRIn) % uncomment for parallel computing
% To reduce the total simulation time, you can execute this loop in
% parallel by using the Parallel Computing Toolbox. Comment out the 'for'
% statement and uncomment the 'parfor' statement. If the Parallel Computing
% Toolbox is not installed, 'parfor' defaults to normal 'for' statement.
% Because parfor-loop iterations are executed in parallel in a
% nondeterministic order, the simulation information displayed for each SNR
% point can be intertwined. To switch off simulation information display,
% set the 'displaySimulationInformation' variable above to false

    % Reset the random number generator so that each SNR point will
    % experience the same noise realization
    rng('default');
    
    % Take full copies of the simulation-level parameter structures so that they are not 
    % PCT broadcast variables when using parfor 
    simLocal = simParameters;
    simLocal.DisplayDiagnostics = false;  % Disable diagnostics in parfor to avoid transparency violations
    waveinfoLocal = waveformInfo;
    
    % Take copies of channel-level parameters to simplify subsequent parameter referencing 
    carrier = simLocal.Carrier;
    pdsch = simLocal.PDSCH;
    pdschextra = simLocal.PDSCHExtension;
    % Create local DL-SCH encoder/decoder objects on the worker to avoid
    % transparency/handle broadcast issues with parfor
    encodeDLSCHLocal = nrDLSCH;
    encodeDLSCHLocal.MultipleHARQProcesses = true;
    encodeDLSCHLocal.TargetCodeRate = pdschextra.TargetCodeRate;

    decodeDLSCHLocal = nrDLSCHDecoder;
    decodeDLSCHLocal.MultipleHARQProcesses = true;
    decodeDLSCHLocal.TargetCodeRate = pdschextra.TargetCodeRate;
    decodeDLSCHLocal.LDPCDecodingAlgorithm = pdschextra.LDPCDecodingAlgorithm;
    decodeDLSCHLocal.MaximumLDPCIterationCount = pdschextra.MaximumLDPCIterationCount;
    decodeDLSCHLocal.reset();        % Reset decoder at the start of each SNR point
    retransInfo = [];
    pathFilters = [];
     
    % Prepare simulation for new SNR point
    SNRdB = simLocal.SNRIn(snrIdx);
    fprintf('\nSimulating transmission scheme 1 (%dx%d) and SCS=%dkHz with %s channel at %gdB SNR for %d 10ms frame(s)\n', ...
        simLocal.NTxAnts,simLocal.NRxAnts,carrier.SubcarrierSpacing, ...
        simLocal.DelayProfile,SNRdB,simLocal.NFrames);
    
    % Specify the fixed order in which we cycle through the HARQ process IDs
    harqSequence = 0:pdschextra.NHARQProcesses-1;

    % Initialize the state of all HARQ processes
    harqEntity = HARQEntity(harqSequence,rvSeq,pdsch.NumCodewords);
    % Initialize a sequential transport-block (TB) counter and a mapping
    % from HARQ process ID -> packet ID so retransmissions can be labeled
    % with the original TB's sequential ID. Index by (procID+1).
    tbCounter = 0;
    tbIDByProc = zeros(pdschextra.NHARQProcesses, pdsch.NumCodewords);
        
    % Create channel locally on the worker and set same properties so each
    % worker has its own instance (avoids transparency violations)
    if contains(simLocal.DelayProfile,'CDL','IgnoreCase',true)
        channelLocal = nrCDLChannel;
        channelLocal = hArrayGeometry(channelLocal,simLocal.NTxAnts,simLocal.NRxAnts);
    else
        channelLocal = nrTDLChannel;
        channelLocal.NumTransmitAntennas = simLocal.NTxAnts;
        channelLocal.NumReceiveAntennas = simLocal.NRxAnts;
    end
    channelLocal.DelayProfile = simLocal.DelayProfile;
    channelLocal.DelaySpread = simLocal.DelaySpread;
    channelLocal.MaximumDopplerShift = simLocal.MaximumDopplerShift;
    channelLocal.SampleRate = waveinfoLocal.SampleRate;
    % Reset the channel instance
    reset(channelLocal);
    
    % Total number of slots in the simulation period
    NSlots = simLocal.NFrames * carrier.SlotsPerFrame;
    
    % Obtain a precoding matrix (wtx) to be used in the transmission of the
    % first transport block (use local channel instance)
    estChannelGridAnts = getInitialChannelEstimate(carrier,simLocal.NTxAnts,channelLocal,simLocal.DataType);
    newWtx = hSVDPrecoders(carrier,pdsch,estChannelGridAnts,pdschextra.PRGBundleSize);
    
    % Timing offset, updated in every slot for perfect synchronization and
    % when the correlation is strong for practical synchronization
    offset = 0;

    % Loop over the entire waveform length
    for nslot = 0:NSlots-1
        
        % Update the carrier slot numbers for new slot
        carrier.NSlot = nslot;
                
        % Calculate the transport block sizes for the transmission in the slot
        [pdschIndices,pdschIndicesInfo] = nrPDSCHIndices(carrier,pdsch);
        trBlkSizes = nrTBS(pdsch.Modulation,pdsch.NumLayers,numel(pdsch.PRBSet),pdschIndicesInfo.NREPerPRB,pdschextra.TargetCodeRate,pdschextra.XOverhead);
        
        % HARQ processing
        for cwIdx = 1:pdsch.NumCodewords
            % If new data for current process and codeword then create a new DL-SCH transport block
            if harqEntity.NewData(cwIdx) 
                trBlk = randi([0 1],trBlkSizes(cwIdx),1);
                setTransportBlock(encodeDLSCH,trBlk,cwIdx-1,harqEntity.HARQProcessID);
                % Assign a new sequential packet ID for this transport block
                tbCounter = tbCounter + 1;
                tbIDByProc(harqEntity.HARQProcessID+1, cwIdx) = tbCounter;
                % If new data because of previous RV sequence time out then flush decoder soft buffer explicitly
                if harqEntity.SequenceTimeout(cwIdx)
                    resetSoftBuffer(decodeDLSCHLocal,cwIdx-1,harqEntity.HARQProcessID);
                end
            end
        end
                                
        % Encode the DL-SCH transport blocks
        codedTrBlocks = encodeDLSCH(pdsch.Modulation,pdsch.NumLayers, ...
            pdschIndicesInfo.G,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
    
        % Get precoding matrix (wtx) calculated in previous slot
        wtx = newWtx;
        
        % Create resource grid for a slot
        pdschGrid = nrResourceGrid(carrier,simLocal.NTxAnts,OutputDataType=simLocal.DataType);
        
        % PDSCH modulation and precoding
        pdschSymbols = nrPDSCH(carrier,pdsch,codedTrBlocks);
        [pdschAntSymbols,pdschAntIndices] = nrPDSCHPrecode(carrier,pdschSymbols,pdschIndices,wtx);
        
        % PDSCH mapping in grid associated with PDSCH transmission period
        pdschGrid(pdschAntIndices) = pdschAntSymbols;
        
        % PDSCH DM-RS precoding and mapping
        dmrsSymbols = nrPDSCHDMRS(carrier,pdsch);
        dmrsIndices = nrPDSCHDMRSIndices(carrier,pdsch);
        [dmrsAntSymbols,dmrsAntIndices] = nrPDSCHPrecode(carrier,dmrsSymbols,dmrsIndices,wtx);
        pdschGrid(dmrsAntIndices) = dmrsAntSymbols;

        % PDSCH PT-RS precoding and mapping
        ptrsSymbols = nrPDSCHPTRS(carrier,pdsch);
        ptrsIndices = nrPDSCHPTRSIndices(carrier,pdsch);
        [ptrsAntSymbols,ptrsAntIndices] = nrPDSCHPrecode(carrier,ptrsSymbols,ptrsIndices,wtx);
        pdschGrid(ptrsAntIndices) = ptrsAntSymbols;

        % OFDM modulation
        txWaveform = nrOFDMModulate(carrier,pdschGrid);

        % Pass data through channel model. Append zeros at the end of the
        % transmitted waveform to flush channel content. These zeros take
        % into account any delay introduced in the channel. This is a mix
        % of multipath delay and implementation delay. This value may 
        % change depending on the sampling rate, delay profile and delay
        % spread
        txWaveform = [txWaveform; zeros(maxChDelay,size(txWaveform,2))]; %#ok<AGROW>
    [rxWaveform,pathGains,sampleTimes] = channelLocal(txWaveform);
        
        % Add AWGN to the received time domain waveform
        % Normalize noise power by the IFFT size used in OFDM modulation,
        % as the OFDM modulator applies this normalization to the
        % transmitted waveform. Also normalize by the number of receive
        % antennas, as the channel model applies this normalization to the
        % received waveform, by default
        SNR = 10^(SNRdB/10);
        N0 = 1/sqrt(simLocal.NRxAnts*double(waveinfoLocal.Nfft)*SNR);
        noise = N0*randn(size(rxWaveform),"like",rxWaveform);
        rxWaveform = rxWaveform + noise;

        if (simLocal.PerfectChannelEstimator)
            % Perfect synchronization. Use information provided by the
            % channel to find the strongest multipath component
            pathFilters = getPathFilters(channelLocal);
            [offset,mag] = nrPerfectTimingEstimate(pathGains,pathFilters);
        else
            % Practical synchronization. Correlate the received waveform
            % with the PDSCH DM-RS to give timing offset estimate 't' and
            % correlation magnitude 'mag'. The function
            % hSkipWeakTimingOffset is used to update the receiver timing
            % offset. If the correlation peak in 'mag' is weak, the current
            % timing estimate 't' is ignored and the previous estimate
            % 'offset' is used
            [t,mag] = nrTimingEstimate(carrier,rxWaveform,dmrsIndices,dmrsSymbols); 
            offset = hSkipWeakTimingOffset(offset,t,mag);
            % Display a warning if the estimated timing offset exceeds the
            % maximum channel delay
            if offset > maxChDelay
                warning(['Estimated timing offset (%d) is greater than the maximum channel delay (%d).' ...
                    ' This will result in a decoding failure. This may be caused by low SNR,' ...
                    ' or not enough DM-RS symbols to synchronize successfully.'],offset,maxChDelay);
            end
        end
        rxWaveform = rxWaveform(1+offset:end,:);

        % Perform OFDM demodulation on the received data to recreate the
        % resource grid, including padding in the event that practical
        % synchronization results in an incomplete slot being demodulated
        rxGrid = nrOFDMDemodulate(carrier,rxWaveform);
        [K,L,R] = size(rxGrid);
        if (L < carrier.SymbolsPerSlot)
            rxGrid = cat(2,rxGrid,zeros(K,carrier.SymbolsPerSlot-L,R));
        end

        if (simLocal.PerfectChannelEstimator)
            % Perfect channel estimation, using the value of the path gains
            % provided by the channel. This channel estimate does not
            % include the effect of transmitter precoding
            estChannelGridAnts = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);

            % Get perfect noise estimate (from the noise realization)
            noiseGrid = nrOFDMDemodulate(carrier,noise(1+offset:end ,:));
            noiseEst = var(noiseGrid(:));

            % Get PDSCH resource elements from the received grid and 
            % channel estimate
            [pdschRx,pdschHest,~,pdschHestIndices] = nrExtractResources(pdschIndices,rxGrid,estChannelGridAnts);

            % Apply precoding to channel estimate
            pdschHest = nrPDSCHPrecode(carrier,pdschHest,pdschHestIndices,permute(wtx,[2 1 3]));
        else
            % Practical channel estimation between the received grid and
            % each transmission layer, using the PDSCH DM-RS for each
            % layer. This channel estimate includes the effect of
            % transmitter precoding
            [estChannelGridPorts,noiseEst] = hSubbandChannelEstimate(carrier,rxGrid,dmrsIndices,dmrsSymbols,pdschextra.PRGBundleSize,'CDMLengths',pdsch.DMRS.CDMLengths); 
            
            % Average noise estimate across PRGs and layers
            noiseEst = mean(noiseEst,'all');

            % Get PDSCH resource elements from the received grid and
            % channel estimate
            [pdschRx,pdschHest] = nrExtractResources(pdschIndices,rxGrid,estChannelGridPorts);

            % Remove precoding from estChannelGridPorts to get channel
            % estimate w.r.t. antennas
            estChannelGridAnts = precodeChannelEstimate(carrier,estChannelGridPorts,conj(wtx));
        end

        % Compute SINR per subcarrier using channel estimates
        gridSize = [carrier.NSizeGrid * 12, carrier.SymbolsPerSlot, pdsch.NumLayers];
        % [k_indices, ~, p_layer] = ind2sub(gridSize, pdschIndices);
        % unique_subcarriers = unique(k_indices);
        % sinr_per_subcarrier = zeros(length(unique_subcarriers), pdsch.NumLayers);
        % for i = 1:length(unique_subcarriers)
        %     subcarrier = unique_subcarriers(i);
        %     for layer = 1:pdsch.NumLayers
        %         mask = (k_indices == subcarrier) & (p_layer == layer);
        %         if any(mask)
        %             h_coeffs = pdschHest(mask, layer);
        %             avg_h_power = mean(abs(h_coeffs).^2);
        %             sinr_per_subcarrier(i, layer) = 10 * log10(avg_h_power / noiseEst);
        %         else
        %             sinr_per_subcarrier(i, layer) = NaN;
        %         end
        %     end
        % end
        % Inputs:
        % gridSize = [carrier.NSizeGrid*12, carrier.SymbolsPerSlot, pdsch.NumLayers];
        % pdschIndices, pdschHest (size: numRE x NumLayers), noiseEst

        [k_idx, ~, p_layer] = ind2sub(gridSize, pdschIndices);   % size: numRE x 1
        nSC = carrier.NSizeGrid * 12;
        nL  = pdsch.NumLayers;

        % Map each RE to the right h_est coefficient
        lin = unique(k_idx);
        v   = abs(pdschHest(lin)).^2;                 % per-RE power

        % If you don’t want average → just keep Psum
        % If you want SINR per RE (not averaged), compute directly:
        sinr_per_re = 10*log10(v / noiseEst);         % vector, same length as pdschIndices

        % Equalization
        [pdschEq,csi] = nrEqualizeMMSE(pdschRx,pdschHest,noiseEst);

        % Common phase error (CPE) compensation
        if ~isempty(ptrsIndices)
            % Initialize temporary grid to store equalized symbols
            tempGrid = nrResourceGrid(carrier,pdsch.NumLayers);

            % Extract PT-RS symbols from received grid and estimated
            % channel grid
            [ptrsRx,ptrsHest,~,~,ptrsHestIndices,ptrsLayerIndices] = nrExtractResources(ptrsIndices,rxGrid,estChannelGridAnts,tempGrid);
            ptrsHest = nrPDSCHPrecode(carrier,ptrsHest,ptrsHestIndices,permute(wtx,[2 1 3]));

            % Equalize PT-RS symbols and map them to tempGrid
            ptrsEq = nrEqualizeMMSE(ptrsRx,ptrsHest,noiseEst);
            tempGrid(ptrsLayerIndices) = ptrsEq;

            % Estimate the residual channel at the PT-RS locations in
            % tempGrid
            cpe = nrChannelEstimate(tempGrid,ptrsIndices,ptrsSymbols);

            % Sum estimates across subcarriers, receive antennas, and
            % layers. Then, get the CPE by taking the angle of the
            % resultant sum
            cpe = angle(sum(cpe,[1 3 4]));

            % Map the equalized PDSCH symbols to tempGrid
            tempGrid(pdschIndices) = pdschEq;

            % Correct CPE in each OFDM symbol within the range of reference
            % PT-RS OFDM symbols
            symLoc = pdschIndicesInfo.PTRSSymbolSet(1)+1:pdschIndicesInfo.PTRSSymbolSet(end)+1;
            tempGrid(:,symLoc,:) = tempGrid(:,symLoc,:).*exp(-1i*cpe(symLoc));

            % Extract PDSCH symbols
            pdschEq = tempGrid(pdschIndices);
        end

        % Decode PDSCH physical channel
        [dlschLLRs,rxSymbols] = nrPDSCHDecode(carrier,pdsch,pdschEq,noiseEst);
        
        % % Display EVM per layer, per slot and per RB
        % if (simLocal.DisplayDiagnostics)
        %     plotLayerEVM(NSlots,nslot,pdsch,size(pdschGrid),pdschIndices,pdschSymbols,pdschEq);
        % end
        
        % Scale LLRs by CSI
        csi = nrLayerDemap(csi); % CSI layer demapping
        for cwIdx = 1:pdsch.NumCodewords
            Qm = length(dlschLLRs{cwIdx})/length(rxSymbols{cwIdx}); % bits per symbol
            csi{cwIdx} = repmat(csi{cwIdx}.',Qm,1);                 % expand by each bit per symbol
            dlschLLRs{cwIdx} = dlschLLRs{cwIdx} .* csi{cwIdx}(:);   % scale by CSI
        end
        
        % Decode the DL-SCH transport channel
        decodeDLSCHLocal.TransportBlockLength = trBlkSizes;
        [decbits,blkerr] = decodeDLSCHLocal(dlschLLRs,pdsch.Modulation,pdsch.NumLayers,harqEntity.RedundancyVersion,harqEntity.HARQProcessID);
        
        % Record retransmission data for each codeword
        for cwIdx = 1:pdsch.NumCodewords
            % Retrieve the packet ID assigned to the HARQ process for this CW
            procID = harqEntity.HARQProcessID; % 0-based
            packetID = tbIDByProc(procID+1, cwIdx);
            % Transmission number (0-based) for this codeword
            retransNum = harqEntity.TransmissionNumber(cwIdx);
            % Compute representative SINR (average over received symbols)
            signalPower = mean(cellfun(@(x) mean(abs(x).^2), rxSymbols));
            noisePower = noiseEst; % noiseEst is a scalar estimate
            sinr_linear = signalPower / noisePower;
            sinr_dB = 10 * log10(sinr_linear);
            retransInfo = [retransInfo, struct('packetID', packetID, 'retransNum', retransNum, 'errorState', blkerr(cwIdx), 'sinr', sinr_dB, 'sinr_per_re', sinr_per_re, 'subcarrier_indices', lin)];
        end
        
        % Store values to calculate throughput
        simThroughput(snrIdx) = simThroughput(snrIdx) + sum(~blkerr .* trBlkSizes);
        maxThroughput(snrIdx) = maxThroughput(snrIdx) + sum(trBlkSizes);
        
        % Update current process with CRC error and advance to next process
        procstatus = updateAndAdvance(harqEntity,blkerr,trBlkSizes,pdschIndicesInfo.G);
        if (simLocal.DisplaySimulationInformation)
            fprintf('\nNSlot=%d/%d, %s',nslot+1,NSlots,procstatus);
        end

        % Get precoding matrix for next slot
        newWtx = hSVDPrecoders(carrier,pdsch,estChannelGridAnts,pdschextra.PRGBundleSize);

     end
    
    % Display the results dynamically in the command window
    if (simLocal.DisplaySimulationInformation)
        fprintf('\n');
    end
    fprintf('\nThroughput(Mbps) for %d frame(s) = %.4f\n',simLocal.NFrames,1e-6*simThroughput(snrIdx)/(simLocal.NFrames*10e-3));
    fprintf('Throughput(%%) for %d frame(s) = %.4f\n',simLocal.NFrames,simThroughput(snrIdx)*100/maxThroughput(snrIdx));
     
    % Store retransmission data for this SNR point
    retransInfoAll{snrIdx} = retransInfo;
end

% Save retransmission data to MAT files after parallel processing
for idx = 1:numel(simParameters.SNRIn)
    retransInfo = retransInfoAll{idx};
    save(fullfile(resultsDir, sprintf('retrans_snr_%d.mat', idx)), 'retransInfo');
end

%% Results
% Display the measured throughput. This is calculated as the percentage of
% the maximum possible throughput of the link given the available resources
% for data transmission.

% Calculate BLER
bler = 1 - simThroughput ./ maxThroughput;

% Save simulation summary to text file
fid = fopen(fullfile(resultsDir, 'simulation_summary.txt'), 'w');
fprintf(fid, 'Simulation Summary\n');
fprintf(fid, '==================\n\n');
fprintf(fid, 'Parameters:\n');
fprintf(fid, 'Delay Profile: %s\n', simParameters.DelayProfile);
fprintf(fid, 'Delay Spread: %.2e s\n', simParameters.DelaySpread);
fprintf(fid, 'Maximum Doppler Shift: %.1f Hz\n', simParameters.MaximumDopplerShift);
fprintf(fid, 'NTxAnts: %d, NRxAnts: %d\n', simParameters.NTxAnts, simParameters.NRxAnts);
fprintf(fid, 'NFrames: %d\n', simParameters.NFrames);
fprintf(fid, 'Total Slots: %d\n', simParameters.NFrames * simParameters.Carrier.SlotsPerFrame);
fprintf(fid, 'Modulation: %s\n', strjoin(string(simParameters.PDSCH.Modulation), '/'));
fprintf(fid, 'HARQ: %s\n\n', simParameters.PDSCHExtension.EnableHARQ);
fprintf(fid, 'Results:\n');
fprintf(fid, 'SNR (dB)\tThroughput (Mbps)\tBLER\n');
for i = 1:length(simParameters.SNRIn)
    fprintf(fid, '%.2f\t\t%.4f\t\t%.2e\n', simParameters.SNRIn(i), 1e-6*simThroughput(i)/(simParameters.NFrames*10e-3), bler(i));
end
fclose(fid);

figure;
plot(simParameters.SNRIn,simThroughput*100./maxThroughput,'o-.')
xlabel('SNR (dB)'); ylabel('Throughput (%)'); grid on;
title(sprintf('%s (%dx%d) / NRB=%d / SCS=%dkHz', ...
              simParameters.DelayProfile,simParameters.NTxAnts,simParameters.NRxAnts, ...
              simParameters.Carrier.NSizeGrid,simParameters.Carrier.SubcarrierSpacing));

% Save the throughput plot
saveas(gcf, fullfile(resultsDir, 'throughput_plot.png'));

% Plot BLER
figure;
semilogy(simParameters.SNRIn, bler, 'o-.')
xlabel('SNR (dB)'); ylabel('BLER'); grid on;
title(sprintf('%s (%dx%d) / NRB=%d / SCS=%dkHz', ...
              simParameters.DelayProfile,simParameters.NTxAnts,simParameters.NRxAnts, ...
              simParameters.Carrier.NSizeGrid,simParameters.Carrier.SubcarrierSpacing));

% Save the BLER plot
saveas(gcf, fullfile(resultsDir, 'bler_plot.png'));

% Bundle key parameters and results into a combined structure for recording
simResults.simParameters = simParameters;
simResults.simThroughput = simThroughput;
simResults.maxThroughput = maxThroughput;
simResults.bler = bler;

%%
% The figure below shows throughput results obtained simulating 10000
% subframes (|NFrames = 1000|, |SNRIn = -18:2:16|).
%
% <<../longRunThroughput.png>>
%

%% Selected Bibliography
% # 3GPP TS 38.211. "NR; Physical channels and modulation."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.212. "NR; Multiplexing and channel coding." 
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.213. "NR; Physical layer procedures for control." 
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TS 38.214. "NR; Physical layer procedures for data."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.
% # 3GPP TR 38.901. "Study on channel model for frequencies from 0.5 to 100
% GHz."
% 3rd Generation Partnership Project; Technical Specification Group Radio
% Access Network.

%% Local Functions

function validateNumLayers(simParameters)
% Validate the number of layers, relative to the antenna geometry

    numlayers = simParameters.PDSCH.NumLayers;
    ntxants = simParameters.NTxAnts;
    nrxants = simParameters.NRxAnts;
    antennaDescription = sprintf('min(NTxAnts,NRxAnts) = min(%d,%d) = %d',ntxants,nrxants,min(ntxants,nrxants));
    if numlayers > min(ntxants,nrxants)
        error('The number of layers (%d) must satisfy NumLayers <= %s', ...
            numlayers,antennaDescription);
    end
    
    % Display a warning if the maximum possible rank of the channel equals
    % the number of layers
    if (numlayers > 2) && (numlayers == min(ntxants,nrxants))
        warning(['The maximum possible rank of the channel, given by %s, is equal to NumLayers (%d).' ...
            ' This may result in a decoding failure under some channel conditions.' ...
            ' Try decreasing the number of layers or increasing the channel rank' ...
            ' (use more transmit or receive antennas).'],antennaDescription,numlayers); %#ok<SPWRN>
    end

end

function estChannelGrid = getInitialChannelEstimate(carrier,nTxAnts,propchannel,dataType)
% Obtain channel estimate before first transmission. This can be used to
% obtain a precoding matrix for the first slot.

    ofdmInfo = nrOFDMInfo(carrier);
    
    chInfo = info(propchannel);
    maxChDelay = chInfo.MaximumChannelDelay;
    
    % Temporary waveform (only needed for the sizes)
    tmpWaveform = zeros((ofdmInfo.SampleRate/1000/carrier.SlotsPerSubframe)+maxChDelay,nTxAnts,dataType);
    
    % Filter through channel    
    [~,pathGains,sampleTimes] = propchannel(tmpWaveform);
    
    % Perfect timing synch    
    pathFilters = getPathFilters(propchannel);
    offset = nrPerfectTimingEstimate(pathGains,pathFilters);
    
    % Perfect channel estimate
    estChannelGrid = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
    
end

function estChannelGrid = precodeChannelEstimate(carrier,estChannelGrid,W)
% Apply precoding matrix W to the last dimension of the channel estimate

    [K,L,R,P] = size(estChannelGrid);
    estChannelGrid = reshape(estChannelGrid,[K*L R P]);
    estChannelGrid = nrPDSCHPrecode(carrier,estChannelGrid,reshape(1:numel(estChannelGrid),[K*L R P]),W);
    estChannelGrid = reshape(estChannelGrid,K,L,R,[]);

end

function plotLayerEVM(NSlots,nslot,pdsch,siz,pdschIndices,pdschSymbols,pdschEq)
% Plot EVM information

    persistent slotEVM;
    persistent rbEVM
    persistent evmPerSlot;
    
    if (nslot==0)
        slotEVM = comm.EVM;
        rbEVM = comm.EVM;
        evmPerSlot = NaN(NSlots,pdsch.NumLayers);
        figure;
    end
    evmPerSlot(nslot+1,:) = slotEVM(pdschSymbols,pdschEq);
    subplot(2,1,1);
    plot(0:(NSlots-1),evmPerSlot,'o-');
    xlabel('Slot number');
    ylabel('EVM (%)');
    legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
    title('EVM per layer per slot');

    subplot(2,1,2);
    [k,~,p] = ind2sub(siz,pdschIndices);
    rbsubs = floor((k-1) / 12);
    NRB = siz(1) / 12;
    evmPerRB = NaN(NRB,pdsch.NumLayers);
    for nu = 1:pdsch.NumLayers
        for rb = unique(rbsubs).'
            this = (rbsubs==rb & p==nu);
            evmPerRB(rb+1,nu) = rbEVM(pdschSymbols(this),pdschEq(this));
        end
    end
    plot(0:(NRB-1),evmPerRB,'x-');
    xlabel('Resource block');
    ylabel('EVM (%)');
    legend("layer " + (1:pdsch.NumLayers),'Location','EastOutside');
    title(['EVM per layer per resource block, slot #' num2str(nslot)]);
    
    drawnow;
    
end