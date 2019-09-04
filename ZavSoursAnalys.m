function ZavSoursAnalys(flNm, proc, isStim, etPrgs, chanMask, chnlGrp, nlxVer)
%ZavSoursAnalys(flNm, proc, isStim, etPrgs, chanMask, chnlGrp, nlxVer)
%read and analysis of source signals (abf-files, Neuralynx, (edf+))
%
%INPUTS
%flNm - name of file (full path)
%proc - type of operation (1 - calculate, 2 - view results)
%isStim - data contain stimulus artefact (if true)
%etPrgs - etalon thresholds for spikes detection
%chanMask - LFP-channels with spikes (without piezo, termo and so on)
%chnlGrp - groups of LFP-channels (not only with spikes)
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%
%OUTPUTS

switch nargin %number of input arguments
    case 0
        disp('file was not asigned')
    case 1
        proc = 1;%calculate
        isStim = false;%no stimulus
        etPrgs = [];%no etalon threshold
        chanMask = [];%no LFP-channel mask
        chnlGrp = {};%no channel groups
        nlxVer = 0;%old version of NLX data (Cheetah 5)
    case 2
        isStim = false;%no stimulus
        etPrgs = [];%no etalon threshold
        chanMask = [];%no LFP-channel mask
        chnlGrp = {};%no channel groups
        nlxVer = 0;%old version of NLX data (Cheetah 5)
    case 3
        etPrgs = [];%no etalon threshold
        chanMask = [];%no LFP-channel mask
        chnlGrp = {};%no channel groups
        nlxVer = 0;%old version of NLX data (Cheetah 5)
    case 4
        chanMask = [];%no LFP-channel mask
        chnlGrp = {};%no channel groups
        nlxVer = 0;%old version of NLX data (Cheetah 5)
    case 5
        chnlGrp = {};%no channel groups
        nlxVer = 0;%old version of NLX data (Cheetah 5)
    case 6
        nlxVer = 0;%old version of NLX data (Cheetah 5)
end

switch proc
    case 1 %calculate
%% === PROCESSING OF THE RECORD === %%
%common procedures for spontaneous and stimulated signals

if all(exist(flNm, 'file') ~= [2, 7])
    [flNm, pth] = uigetfile({'*.abf'; '*.ncs; *.nev; *.nse'}, 'Select file', flNm);%open dialog
    flNm = [pth, flNm];%full pathname
    if (all(flNm == 0) || all(pth == 0))
        return;
    end
end

bufStr = flNm(end - 3:end);%extension
if (exist(flNm, 'dir') == 7) %directory (Neuralynx)
    if (flNm(end) ~= '\') %no slash
        flNm(end + 1) = '\';
    end
    svFlNm = [flNm(1:end - 1), '.mat'];%pathname of file for calculated values
else %((exist(flNm, 'file') == 2) && (strcmp(bufStr, '.abf') || strcmp(bufStr, '.edf') || strcmp(bufStr, '.daq')))%file ("abf" or "edf")
    svFlNm = [flNm(1:end - 3), 'mat'];%pathname of file for calculated values
end
zavp = struct('file', [], 'siS', [], 'dwnSmplFrq', [], 'stimCh', [], 'realStim', [], 'rarStep', []);%structure with parameters
zavp.file = flNm;%file path and name

largestF = 1;%number of file (in ncsFiles structure) with largest size
if ((exist(flNm, 'dir') == 7) && ~strcmp(bufStr, '.abf') && ~strcmp(bufStr, '.edf') && ~strcmp(bufStr, '.daq')) %Neuralynx-files
    dirCnt = dir(flNm);
    for t = 1:length(dirCnt)
        if (~dirCnt(t).isdir && ~isempty(strfind(dirCnt(t).name, '.ncs')))
            if (dirCnt(t).bytes > dirCnt(largestF).bytes)
                largestF = t;%number of largest neuralynx-file
            end
        end
    end
    largestF = str2double(dirCnt(largestF).name(4:(end - 4)));%number of channel in Neuralynx-system
%else%for other types of file
%    largestF = 1;%first channel
end
disp(flNm)

[data, synchro, hd] = ZavLoadData(flNm, [], largestF, 1, [], [], nlxVer);%load header of file
[zavp.stimCh, zavp.realStim] = StimMoments(flNm, hd, isStim, synchro);%define moments of stimulus
%data = ZavLoadData(flNm, hd, largestF, 1, 0, 'e', nlxVer);%one channel and one sweep
recLen = size(data, 1);%memory preallocation (length of single-sweep largest recordation)

zavp.siS = hd.si * 1e-6;%sampling interval (seconds)
zavp.dwnSmplFrq = 1000;%downsampled-discretization frequency (Hz, for signals to be treated)

if ((hd.nOperationMode == 3) && (hd.lActualEpisodes == 0))%data were acquired in gap-free mode
    hd.lActualEpisodes = 1;%single sweep for gap-free mode (abf)
end

fs = 1e6 / hd.si;%sampling frequency (Hz)
if (fs > zavp.dwnSmplFrq) %downsampling
    bfrr = resample(data, zavp.dwnSmplFrq, fs);%LFP resampled down to zavp.dwnSmplFrq Hz
else %upsampling
    tmOrig = (0:hd.si:((size(data, 1) - 1) * hd.si))';%original time
    tmResm = (0:1e3:((size(data, 1) - 1) * hd.si))';%resempled (downsampled) time
    bfrr = interp1(tmOrig, data, tmResm);%resampled lfp
end
%interval corresponding to 1(one) millisecond ((1e6 / zavp.dwnSmplFrq) microseconds)
zavp.rarStep = length(data) / length(bfrr);%step of rarefy to get 1 ms descretezation (samples)

lfp = zeros(size(bfrr, 1), hd.nADCNumChannels, hd.lActualEpisodes);%downsampled lfp %hd.sweepLengthInPts - length of one sweep on a channel
spks(1:hd.nADCNumChannels, 1:hd.lActualEpisodes) = struct('tStamp', [], 'ampl', [], 'shape', []);%moments of spikes appearance (ms)
lfpVar = Inf(hd.nADCNumChannels, hd.lActualEpisodes);%minimal std of lfp (lowest noise level in the record)
spkTime = ((-1e-3:zavp.siS:1e-3) / zavp.siS)';%samples of spike time course (-1ms to 1ms from peak)


% if (isfield(hd, 'dspDelay_mks') && any(diff(hd.dspDelay_mks) ~= 0)) %different delays
%     disp('=============================================================================')
%     disp(flNm)
%     disp('We have to compensate difference in dspDelay!')
%     disp('=============================================================================')
%     return;
% end
    %     %time shift!
    %     if (hd.dspDelay_mks(ch) == 0)
    %         data = [zeros(16, 1); data(1:(end - 16))];
    %     end

    
tstp = 5e5;%time step (samples)
if isempty(chnlGrp) %no channel groups indicated
    shnkChN = 16;%32;%default number of channels in a shank
    chnlGrp = cell(1, ceil(hd.nADCNumChannels / shnkChN));
    ch = 1;
    for t = 1:length(chnlGrp) %run over N-channels groups
        chnlGrp(t) = {ch:min((ch + (shnkChN - 1)), hd.nADCNumChannels)};%groups of channels (shanks)
        ch = chnlGrp{t}(end) + 1;%first channel of next shank
    end
end
for t = 1:length(chnlGrp)
    disp(['chgrp', num2str(t), ': ', num2str(chnlGrp{t}(1)), '-', num2str(chnlGrp{t}(end))])
end
for cg = 1:length(chnlGrp) %run over channel groups (shanks)
    for sw = 1:hd.lActualEpisodes %run over sweeps (segments, episodes)
        dataFltCmn = zeros(recLen, length(chnlGrp{cg}), 'single');%filtered traces from all channels of current shank
        mednCh = zeros(recLen, 1, 'single');%spatial median (over channels)
        dsprCh = zeros(recLen, 1, 'single');%spatial dispertion (over channels)
        
        %(1) common reference calculation
        for nk = 1:length(chnlGrp{cg}) %run over channel in the group
            ch = chnlGrp{cg}(nk);%real number of channel
%             disp(['read channel ', num2str(ch), ' of ', num2str(hd.nADCNumChannels)])
            data = ZavLoadData(flNm, hd, ch, sw, [], [], nlxVer);%read one channel all sweeps (no nlx-spikes)
            data = squeeze(data);%delete excess dimensions
            if (size(data, 1) == recLen) %full length recordation
                [pntsL, pntsR] = ZavFindFlats(data);
                jj = (diff([pntsL, pntsR], 1, 2) <= zavp.rarStep);%number of short segment (<1ms)
                pntsL(jj) = []; pntsR(jj) = [];%1ms and longer silent segments only
                if (zavp.rarStep >= 1) %downsampling
                    bfrr = resample(data, zavp.dwnSmplFrq, fs);%LFP downsampled to zavp.dwnSmplFrq Hz
                    lfp(1:numel(bfrr), ch, sw) = bfrr;%downsampled LFP
                else %upsampling
                    lfp(:, ch, sw) = interp1(tmOrig, data, tmResm);%upsampled LFP
                end
                if (zavp.siS < 1e-3) %fast discretisation necessary
                    dataFlt = ZavFilter(data, fs, 'high', 300, 2);%highpass filtration
                    if ((fs / 2.001) > 5000) %fast discretisation
                        dataFlt = ZavFilter(dataFlt, fs, 'low', 5000, 2);%lowpass filtration
                    end

                    %differential-filter (instead of conventional)
                    %dataFlt = [diff(data); 0];%differential-filter
                else %slow discretisation
                    dataFlt = ZavFilter(data, fs, 'high', round(fs / 2.5), 2);%highpass filtration
                end
                lfpVar(ch, sw) = CalcMinVar(dataFlt, zavp.siS, data, 1);%lfp(:, ch, sw));%minimal dispers on a trace
                dataFltCmn(:, nk) = single(dataFlt);%write current channel
                for t = 1:length(pntsL) %run over flat segments
                    jj = pntsL(t):pntsR(t);%points of flat segment
                    dataFltCmn(jj, nk) = NaN;%replace flat segments with NaN
                end
            else
%                 disp('skipped channel')
                dataFltCmn(:, nk) = NaN;%write current channel
            end
        end
        if (length(chnlGrp{cg}) > 2) %many independent channels recorded
            for t = 1:tstp:recLen %run over time
                jj = t:min((t + tstp - 1), recLen);
                dsprCh(jj) = nanstd(dataFltCmn(jj, :), [], 2);%spatial dispertion
                mednCh(jj) = nanmedian(dataFltCmn(jj, :), 2);%spatial median
            end
        end
%         bins = 0:0.1:30;%binning for distribution
%         [~, t] = max(histc(dsprCh, bins));%peak of distribution
%         z = bins(t);%most frequent value
%         aa = dsprCh((dsprCh <= z) & (dsprCh > 0));
%         aa = [aa; (2 * z) - aa];
%         dspPrg = z + 4 * std(aa);%threshold for deviation of dispertion

        %(2) treat filtered trace (LFP deviation, spikes detection uning common reference)
        for nk = 1:length(chnlGrp{cg}) %run over channels in the group
            ch = chnlGrp{cg}(nk);%real number of channel
            if (ismember(ch, chanMask) && (zavp.siS < 1e-3)) %only channels contaning EEG-signal (fast discretisation necessary)
%                 disp(['spikes of channel ', num2str(ch)])
                if ~isempty(etPrgs) %alternative threshold assigned by user
                    etTrshld = etPrgs(ch);%etalon thresholds for spikes detection
                end
                if isempty(etPrgs) %no etalon thresholds provided
                    etTrshld = -3 * lfpVar(ch, sw);%newly calculated threshold
                end

                %common reference
                dataFlt = double(dataFltCmn(:, nk) - mednCh);%subtract common reference

                [spks(ch, sw).tStamp, spks(ch, sw).ampl] = ZavFindSpikes(dataFlt, fs, etTrshld, false, false);%moments of peaks of spikes appearance (samples (from 1))
                spks(ch, sw).tStamp = (spks(ch, sw).tStamp - 1) / zavp.rarStep;%moments of peaks of spikes appearance (ms from record beginning (from 0))
                %[~, ~, ~, spkNlx] = ZavLoadData(flNm, hd, ch, sw, [], '[], nlxVer);%read one channel all sweeps (plus nlx-spikes)
                %spks(ch, sw).tStamp = spkNlx.tStamp(:) / 1e3;%Neuralynx spikes: mks to ms (colon!)
            end
%             %get shapes of spikes
%             if ~isempty(spks(ch, sw).tStamp) %spikes exist
%                 spks(ch, sw).tStamp(round((spks(ch, sw).tStamp * zavp.rarStep) + spkTime(1)) < 1) = [];%delete spikes on edge
%                 spks(ch, sw).tStamp(round((spks(ch, sw).tStamp * zavp.rarStep) + spkTime(end)) > size(data, 1)) = [];%delete spikes on edge
%                 dataFlt = ZavFilter(data, 1 / zavp.siS, 'high', 300, 2);%filtration (mode 2 - cascade (double) filtration)
%                 jj = repmat(spks(ch, sw).tStamp' * zavp.rarStep, numel(spkTime), 1);%complite indices for spikes extraction
%                 jj = round(jj + repmat(spkTime, 1, numel(spks(ch, sw).tStamp)));%complite indices for spikes extraction
%                 spks(ch, sw).shape = dataFlt(jj);%spikes shape
%             end
        end
    end
end
clear dataFltCmn dsprCh mednCh data
%set adequate precision of LFP
lfp = single(lfp);%conver to single precision
for ch = 1:hd.nADCNumChannels %run over channels
    for sw = 1:hd.lActualEpisodes %run over sweeps (segments, episodes)
        spks(ch, sw).ampl = single(spks(ch, sw).ampl);%conver to single precision
        spks(ch, sw).shape = single(spks(ch, sw).shape);%conver to single precision
   end
end
save(svFlNm, 'hd', 'zavp', 'lfp', 'spks', 'lfpVar', 'chnlGrp', '-v7.3')%saving calculated values

% %false spikes separation
% for sw = 1:hd.lActualEpisodes %run over sweeps (segments, episodes)
%     for cg = 1:length(chnlGrp) %run over channel groups (shanks)
%         spks(:, sw) = ZavCheckSpikes(spks(:, sw), lfpVar, chnlGrp{cg}, 12, 0.2, size(lfp, 1));%discard unreal spikes
%     end
% end

% %set adequate precision spikes
% for ch = 1:hd.nADCNumChannels %run over channels
%     for sw = 1:hd.lActualEpisodes %run over sweeps (segments, episodes)
%         spks(ch, sw).ampl = single(spks(ch, sw).ampl);%conver to single precision
%         spks(ch, sw).shape = single(spks(ch, sw).shape);%conver to single precision
%    end
% end
% save(svFlNm, 'spks', '-append')%saving calculated values


%%% === VIEW RESULTS === %%%
    case 2 %view results
        
load(flNm)%load results of data treatment


%% === VIEW SIGNALS === %%
lfp = double(lfp);
for sw = 1:hd.lActualEpisodes
    for ch = 1:hd.nADCNumChannels
        jj = spks(ch, sw).ampl <= (-5 * lfpVar(ch));
        spks(ch, sw).tStamp = spks(ch, sw).tStamp(jj);
        spks(ch, sw).ampl = spks(ch, sw).ampl(jj);
    end
end

%segms(:, 1) - position of synchro-event (ms from record begin (from 0))
%segms(:, 2) - number of channel where syncro-event was detected
%segms(:, 3) - number of sweeps where syncro-event was detected
%segms(:, 4) - number of stimulus (within sweep, in matrix zavp.realStim) or number of trough in brst matrix

segms = zeros(numel(vertcat(zavp.realStim(:).r)), 4);%preallocation of memory for stimuli moments
z = 1;%through counter of stimuli
for sw = 1:hd.lActualEpisodes
    for t = 1:numel(zavp.realStim(sw).r) %run over synchro-events
        segms(z, 1) = zavp.realStim(sw).r(t) / zavp.rarStep;%position of synchro-event (ms from record begin (from 0))
        %segms(z, 1) = segms(z, 1) + sepOnsetPeak2(7, sw).r(z, 2); disp('. shifted start point .') %shift to SEP onset
        segms(z, 2) = zavp.stimCh;%number of channel where syncro-event was detected
        segms(z, 3) = sw;%number of sweep where syncro-event was detected
        segms(z, 4) = t;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
        z = z + 1;%through counter of stimuli
    end
end
segms(z:end, :) = [];%delete excess

if (isempty(segms) || isequal(hd.fFileSignature, 'eeg EDF+'))%isempty(zavp.realStim(1).r)
    tW = 5e3;%time window
    tm = (0:tW:size(lfp, 1))';
    segms = [tm, repmat([NaN, 1], length(tm), 1), (1:length(tm))'];
%     if exist([zavp.file(1:end - 1), '_brst.mat'], 'file') %bursts
%         load([zavp.file(1:end - 1), '_brst.mat'])
%         tm = vertcat(brst(:).t); ch = brst(1).ch(1);
%         segms = [tm(:, 1), repmat([ch, 1], length(tm), 1), (1:length(tm))'];
%     elseif exist([zavp.file(1:end - 4), '_mbs.mat'], 'file') %mechano-movement detected (mbs)
%         load([zavp.file(1:end - 4), '_mbs.mat'])
%         segms = [mbs(:, 1), repmat([mdCh(1), 1], size(mbs, 1), 1), (1:size(mbs, 1))'];
%     end
end
% lenR = [0, Inf];%movement length range (larger than lenR(1) and shorter than lenR(2))    

% set(gcf, 'DefaultTextFontSize', 16, 'DefaultTextFontName', 'Arial', 'DefaultAxesFontSize', 16, 'DefaultAxesFontName', 'Arial');
% set(gcf, 'Units', 'normalized', 'Position', [0.003, 0.046, 0.45391, 0.3916]) %compact
% % set(gcf, 'Units', 'normalized', 'Position', [0.003, 0.046, 1, 0.3916]) %horizontal
% % set(gcf, 'Units', 'normalized', 'Position', [0.003, 0.046, 0.45391, 0.88379]) %vertical
%
% % set(gcf, 'Units', 'normalized', 'Position', [1.003, 0.003, 1, 0.3916]) %second monitor, horizontal
% % set(gcf, 'Units', 'normalized', 'Position', [1.003, 0.003, 0.45391, 0.88379]) %second monitor, vertical

hh = floor(hd.recTime / 3600);%hours
mm = floor((hd.recTime / 60) - hh * 60);%minutes
ss = hd.recTime - hh * 3600 - mm * 60;%seconds
disp(['from    ', num2str(hh(1)), ':', num2str(mm(1)), ':', num2str(ss(1))])
disp(['to      ',num2str(hh(2)), ':', num2str(mm(2)), ':', num2str(ss(2))])

%power spectrum parameters (chronux)
clear prm
if ~isfield(zavp, 'dwnSmplFrq') %field exist
    zavp.dwnSmplFrq = 1e3;%discretization frequency for signals to be treated (Hz)
end
prm.Fs = zavp.dwnSmplFrq;%sampling frequency for resampled data
prm.fpass = [25, 70];%frequency range (predominantly for sensory response)
prm.tapers = [2, 3];%tapers
prm.pad = 1;%padding factor
prm.err = [2, 0.05];
prm.trialave = 0;%channel by channel, sweep by sweep (don't average when 0)
%oscillations frequency ranges% [1, 4] - delta; [4, 8] - theta;
%                               [8, 13] - alpha; [13, 30] - beta; [30, 80] - gamma


%% === LFP, CSD, MUA === %%
nlxVer = 1;%0 - old format of Cheetah log-file (1 - new format)
rawData = 0;%load raw data if true

pltCSD = 0;%plot CSD if true
pltMUA = 0;%plot MUA if true

dCoeff(1) = -1.0; dCoeff(2) = 10000;%max(max(abs(whtLFP)));%max(max(abs(lfp(:, rCh((n + 1):end)))));%characteristic sizes (amplitude coefficient)
caRng = 500;%range of color axis for CSD %max(max(abs(lfpCSD)));

segmEdge = [0, 200e3] + 0;%left and right shifts from synchro-point (ms)
rCh = 1;%:16;%channels to be read and draw (1:hd.nADCNumChannels)
%EEG: [30, 1, 5, 9, 13, 17, 21, 25]; %[30, 4, 8, 16, 20, 24, 27];
%poly3: (1:10); %(11:22); %(23:32);
tecCh = 0;%number of special channels

sn = 12;%:size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));
if isempty(sn), return; end

if rawData %raw data requested
    prm.Fs = 1 / zavp.siS;%sampling frequency for raw data
else %resampled data requested
    prm.Fs = zavp.dwnSmplFrq;%sampling frequency for rare data
end
tm = segmEdge(1):(1e3 / prm.Fs):segmEdge(2);%time matrix for raw or downsampled data (ms)

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, rCh, rawData, nlxVer);%lfp phased with respect to synchro-events
% whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to synchro-events

whtLFP = whtLFP - repmat(median(whtLFP((tm >= -0) & (tm <= 50), :, :), 1), size(whtLFP, 1), 1);%remove DC
% whtLFP = whtLFP - repmat(median(whtLFP((tm >= -100) & (tm <= -85), :, :), 1), size(whtLFP, 1), 1);%remove DC
% whtLFP = whtLFP - repmat(median(whtLFP(1:1e2, :, :), 1), size(whtLFP, 1), 1);%remove DC

% wLFPcomm = mean(whtLFP(:, 1:(end - tecCh), :), 2);%common reference
wLFPcomm = mean(ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, 1:16, rawData, nlxVer), 2);%common reference
% wLFPcomm = mean(mean(ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, 1:16, rawData, nlxVer), 2), 3);%common reference

%common reference
whtLFP = whtLFP - repmat(wLFPcomm, 1, size(whtLFP, 2), 1);%common reference
    
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)


%special channels
for ch = numel(rCh) + (0:-1:(1 - tecCh))
    whtLFP(:, ch) = (whtLFP(:, ch) - mean(whtLFP(:, ch))) * 0.5;%scalling of signal (piezo)
end

%filtration
for ch = 1:(size(whtLFP, 2) - tecCh) %run over LFP-channels
    %0
%     jj = 1:1e3:(2 * 60e3);
%     p1 = polyfit(tm(jj), whtLFP(jj, ch), 1);
%     whtLFP(:, ch) = whtLFP(:, ch) - (p1(1) * tm + p1(2));
    %1
%     whtLFP(:, ch) = locdetrend(whtLFP(:, ch), 1e3, [300, 10]);
    %2
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'high', 1, 2);%
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'low', 5000, 2);%
%     whtLFP(:, ch) = locdetrend(whtLFP(:, ch), prm.Fs, [0.002, 0.001]);%for raw data (spikes)
    %3
%     whtLFP(:, ch) = ZavSpctSubtr(whtLFP(:, ch), []);%denoise by spectral subtraction
    
    %4 (RC-filter)
%     whtLFP(:, ch) = ZavRCfilt(whtLFP(:, ch), 2, prm.Fs, 'high');
end

% if (diff(segmEdge) < 5e3)
%     dCoeff(4) = 1; bufStr = 'ms';%time coefficient
% elseif ((diff(segmEdge) >= 5e3) && (diff(segmEdge) < 5e4))
    dCoeff(4) = 1e3; bufStr = 's';%time coefficient
% else
%     dCoeff(4) = 60e3; bufStr = 'min';%time coefficient
% end

clf, hold on
if (pltCSD && (size(whtLFP, 2) > 5)) %plot CSD for many channels
    ajj = whtLFP(:, 1:(end - tecCh));
%     if (size(whtLFP, 1) > 1e3)
% %         ajj = ZavFilter(whtLFP, 1e3, 'high', 5, 2);%remove DC
%         ajj = ZavRCfilt(whtLFP, 1, 1e3, 'high');%remove DC
%     else
% %         ajj = whtLFP - repmat(mean(whtLFP((tm >= 1) & (tm <= 3), :), 1), numel(tm), 1);%remove DC
%         ajj = whtLFP - repmat(mean(whtLFP(1:min(10, size(whtLFP, 1)), :), 1), numel(tm), 1);%remove DC
%     end
%     ajj = ajj ./ repmat(mean(lfpVar(rCh, :), 2)' / min(mean(lfpVar(rCh, :), 2)), numel(tm), 1);%amplitude correction on impedance
% %     ajj = filtfilt([0.607, 1, 0.607], 2.213, ajj')';%spatial filtration    

    lfpCSD = -diff(ajj, 2, 2);%CSD (current source density)
    k = 5;%round(diff(segmEdge) / 100) + 1;%interval of smoothing
    for t = 1:size(lfpCSD, 2) %run over (channels) depth
        lfpCSD(:, t) = smooth(lfpCSD(:, t), k);%smoothing along time
    end
    %lfpCSD = filtfilt([0.607, 1, 0.607], 2.213, lfpCSD')';%3-channel spatial smoothing as in elephy
    
    %interpolation along channels
    ajj = lfpCSD;
    spPrc = 0.1;%1;%spatial precision
    lfpCSD = zeros(size(ajj, 1), numel(1:spPrc:((numel(rCh) - tecCh) - 2)));
    for t = 1:size(lfpCSD, 1) %run over time
        lfpCSD(t, :) = interp1(1:((numel(rCh) - tecCh) - 2), ajj(t, :), 1:spPrc:((numel(rCh) - tecCh) - 2));%interpolation along channels
        %lfpCSD(t, :) = spline(1:((numel(rCh) - spchn)- 2), ajj(t, :), 1:spprc:((numel(rCh) - spchn) - 2));%spline interpolation along channels
    end
    
    imagesc(tm / dCoeff(4), 2:((numel(rCh) - tecCh) - 1), lfpCSD')
    set(gca, 'Color', 0.5 * ones(1, 3))%%gray rectangle background %imagesc(tm([1, end]) / dCoeff(4), [1, numel(rCh)], -caRng * ones(size(lfpCSD, 2), 2))
    colormap('default');
    hclrb = colorbar;%colobar handle
    set(get(hclrb, 'YLabel'), 'String', 'sink - source')

    caxis(caRng * [-1, 1])%colorbar range
    %caxis(max(max(abs(lfpCSD))) * [-1, 1]);%colorbar range
end

%= plot LFP =%
% if (size(whtLFP, 1) > 3e4)
%     ii = round(linspace(1, size(whtLFP, 1), 3e4));%rar index
% else
    ii = 1:size(whtLFP, 1);%full rastr
% end
dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient
plot(tm(ii) / dCoeff(4), (whtLFP(ii, :) / dCoeff(3)) + repmat(1:numel(rCh), numel(ii), 1), 'k', 'LineWidth', 1)
%= end of plot LFP =%

%= plot MUA =%
if ((numel(sn) == 1) && pltMUA) %single sweep (segment)
    jj = segms(sn, 1) + segmEdge;%absolute time of the segment boundary ( ms, e.i. from begin of file)
    sw = segms(sn, 3);%sweep
    for ch = 1:numel(rCh) %run over wanted channels
        spkDens = spks(rCh(ch), sw).tStamp((spks(rCh(ch), sw).tStamp > jj(1)) & (spks(rCh(ch), sw).tStamp < jj(end))) - segms(sn, 1);%unites in the segment
        spkCntDens = NaN(3 * length(spkDens), 2);%continuous trace of spikes
        spkCntDens(1:3:end, 1) = spkDens;%time point of lower level
        spkCntDens(2:3:end, 1) = spkDens;%time point of upper level
        spkCntDens(1:3:end, 2) = ch - 0.4;%lower level
        spkCntDens(2:3:end, 2) = ch - 0.2;%upper level
        % %plot([spkDens'; spkDens'] / dCoeff(4), ch + repmat([0; 0.2] - 0.4, 1, numel(spkDens)), 'r', 'LineWidth', 1.5)
        plot(spkCntDens(:, 1) / dCoeff(4), spkCntDens(:, 2), 'r', 'LineWidth', 1.5) %fast drawing of spikes
    end
end
%= end of plot MUA =%

set(gca, 'YDir', 'reverse', 'XLim', segmEdge / dCoeff(4), 'YLim', [-0.05, numel(rCh) + 1.05])
set(get(gca, 'XLabel'), 'String', ['time, ', bufStr])

%disp([num2str((segms(sn([1, end]), 1)' + segmEdge) / 1e3), ' s'])
if (size(hd.recChNames, 2) > 1) %alternames was signed
    chnlStr = cell(length(rCh), 1);
    for t = 1:length(rCh) %run over channel
        chnlStr{t} = [num2str(rCh(t)), '(', hd.recChNames{rCh(t), 2}, ')'];%alter name of channels
    end
else
    chnlStr = num2str(rCh(:));%numbered channels
end
set(gca, 'YTick', 1:numel(rCh), 'YTickLabel', chnlStr)

set(get(gca, 'YLabel'), 'String', 'channels', 'Units', 'normalized', 'Position', [-0.06, 0.6, 1])
plot([0, 0], get(gca, 'YLim'), 'm', 'LineWidth', 2)%nul-time line

[~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
annotation('textarrow', 'Position', [0.05, 0.15, 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    'String', [num2str(round(dCoeff(2) * 10) / 10), ' \muV'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];%many segments
else
    bufStr = [bufStr, num2str(sn)];%single segment
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if (isempty(p2) || (p2 < p1)), p2 = length(zavp.file) - 1; end
fN = [zavp.file(p1:p2), bufStr, '_LFP'];
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.eps'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title([fN, ' ', num2str((segms(sn(1), 1) + segmEdge(1)) / 1e3), '-', num2str((segms(sn(end), 1) + segmEdge(2)) / 1e3), ' s'])



%% === MUA density === %%
segmEdge = [-50e0, 50e0] + 0;%left and right shifts from synchro-point (ms)
dCoeff(1) = -1; dCoeff(2) = 5;%50;%max(max(spkDens));%characteristic sizes
caRng = 5;%range of MUA colormap %max(max(spkDens));%
clr = 'k';%'w' - show colored spike density; any else - plot timecourses of spike density 

rCh = [2, 13];%1:16;%channels to be read and draw (1:hd.nADCNumChannels)
%EEG: [30, 1, 5, 9, 13, 17, 21, 25]; %[30, 4, 8, 16, 20, 24, 27];
%poly3: (1:10); %(11:22); %(23:32);

sn = 20;%:size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));
if isempty(sn), disp('no one segment choosed'); return; end

%tm = segmEdge(1):segmEdge(2);%time matrix for raw or downsampled data (ms)
prm.Fs = 1 / (zavp.siS * zavp.rarStep);%sampling frequency for rare data

binStep = [100, 5];%[5, 5];%[bin size (ms), precision (must be integer)]
[~, basbins] = ZavOverHist(spks(rCh(1), sw).tStamp, binStep, segmEdge);%fast overlaped histogram (histc)
spkDens = zeros(length(basbins), numel(rCh), numel(sn));%memory preallocation for density of spikes in time (histogram)

for z = 1:numel(sn) %run over requested segments
    jj = segms(sn(z), 1) + segmEdge;%time segment boundary (ms from begin of file)
    sw = segms(sn(z), 3);%sweep
    for ch = 1:numel(rCh) %run over channels
        spkDens(:, ch, z) = ZavOverHist(spks(rCh(ch), sw).tStamp, binStep, jj);%fast overlaped histogram (histc)
    end
end
spkDens = mean(spkDens, 3) / binStep(1);%mean spikes density

% %delete artificial spikes near stimuli
% ii = (basbins > -2) & (basbins < 1); spkDens(ii, :) = 0;%erase artefacts

% if (diff(segmEdge) < 5e3)
    dCoeff(4) = 1; bufStr = 'ms';%time coefficient
% elseif ((diff(segmEdge) >= 5e3) && (diff(segmEdge) < 5e4))
%     dCoeff(4) = 1e3; bufStr = 's';%time coefficient
% else
%     dCoeff(4) = 60e3; bufStr = 'min';%time coefficient
% end

clf, hold on
if (numel(rCh) > 1) %show many channels
%     k = 3;%round(diff(segmEdge) / 100) + 1;%interval of smoothing
%     for t = 1:size(spkDens, 2)
%         spkDens(:, t) = smooth(spkDens(:, t), k);%smooth along time
%     end
    
    if (clr == 'w') %show colored spike density
        ajj = spkDens;%non-smoothed colormap
        %interpolation along channels (smoothing)
%         ajj = zeros(size(spkDens, 1), numel(1:0.25:numel(rCh)));%wiil contain density interpolated along channels
%         for t = 1:size(spkDens, 1)
%             ajj(t, :) = interp1(1:numel(rCh), spkDens(t, :), 1:0.25:numel(rCh));%interpolation along channels
%         end
        imagesc(basbins / dCoeff(4), 1:numel(rCh), ajj')%interpolated
        
%         imagesc(basbins / dCoeff(4), 1:numel(rCh), spkDens')
        %clrMap = colormap('gray'); colormap(clrMap(end:-1:1, :))
        shading flat; colormap('default');
        hclrb = colorbar;
        set(get(hclrb, 'YLabel'), 'String', ['spikes/', num2str(binStep(1)), 'ms'])
        caxis([0, caRng])%colorbar range
    end
    set(gca, 'YDir', 'reverse')%orientation

    %plot spike density curves
    dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient
    
%     %average on channels
%     ii = (bins < 5);
%     spkDens = spkDens - repmat(mean(spkDens(ii, :), 1), size(spkDens, 1), 1);%base minus
%     ch = {1:2, 5:6, 8, 10:11, 13:14};
%     spkDens1 = [];
%     for t = 1:numel(ch), spkDens1 = [spkDens1, mean(spkDens(:, ch{t}), 2)]; end
%     rCh = 1:numel(ch);
%     spkDens = spkDens1;
    
    plot(basbins / dCoeff(4), (spkDens / dCoeff(3)) + repmat(1:numel(rCh), numel(basbins), 1), clr, 'LineWidth', 1.5)
    set(gca, 'YTick', 1:numel(rCh), 'YTickLabel', num2str(rCh'), 'YLim', [0, numel(rCh) + 1])
    [~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
    annotation('textarrow', 'Position', [0.05, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
        'String', [num2str(round(dCoeff(2) * 10) / 10), 'spikes/', num2str(binStep(1)), 'ms'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');
    set(get(gca, 'YLabel'), 'String', 'channels')
elseif (numel(rCh) == 1) %show single channel
    plot(basbins / dCoeff(4), spkDens, clr, 'LineWidth', 2)
    set(get(gca, 'YLabel'), 'String', ['spiks/', num2str(binStep(1)), 'ms'])
end
plot([0, 0], get(gca, 'YLim'), 'm', 'LineWidth', 2)%nul-time line
       
set(get(gca, 'XLabel'), 'String', ['time, ', bufStr])
set(gca, 'XLim', segmEdge / dCoeff(4))

%disp([num2str((segms(sn([1, end]), 1)' + segmEdge) / 1e3), ' s'])
if (size(hd.recChNames, 2) > 1) %alternames was signed
    chnlStr = hd.recChNames(rCh, 2);%alter name of channels
else
    chnlStr = num2str(rCh(:));%numbered channels
end
set(gca, 'YTick', 1:numel(rCh), 'YTickLabel', chnlStr)

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];%many segments
else
    bufStr = [bufStr, num2str(sn)];%single segment
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if (isempty(p2) || (p2 < p1)), p2 = length(zavp.file) - 1; end
fN = [zavp.file(p1:p2), bufStr, '_MUA'];
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.eps'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title([fN, ' ', num2str((segms(sn(1), 1) + segmEdge(1)) / 1e3), '-', num2str((segms(sn(end), 1) + segmEdge(2)) / 1e3), ' s'])


%% === LFPfromMUA, CSDfromMUA, MUA === %%
rawData = 0;%raw data (1 - load raw)
pltCSD = 1;%plot CSD (1 - plot)
dCoeff(1) = -1.0; dCoeff(2) = 5;%max(max(abs(whtLFP)));%max(max(abs(lfp(:, rCh((n + 1):end)))));%characteristic sizes (amplitude coefficient)
caRng = 10;%range of color axis for CSD %max(max(abs(lfpCSD)));

segmEdge = [0, 50];%left and right shifts from synchro-point (ms)
rCh = 1:16;%channels to be read and draw (1:hd.nADCNumChannels)
%EEG: [30, 1, 5, 9, 13, 17, 21, 25]; %[30, 4, 8, 16, 20, 24, 27];
%poly3: (1:10); %(11:22); %(23:32);

sn = 1:150;%size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));
if isempty(sn), return; end

if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k))' / k;%time matrix for raw or downsampled data (ms)

whtLFP = zeros(10 * size(tm, 1), length(rCh), length(sn));
for z = 1:numel(sn) %run over requested segments
    jj = segms(sn(z), 1) + segmEdge;%time segment boundary (ms from begin of file)
    sw = segms(sn(z), 3);%sweep
    for ch = 1:numel(rCh) %run over channels
        ii = (spks(rCh(ch), sw).tStamp >= jj(1)) & (spks(rCh(ch), sw).tStamp <= jj(end));%spikes within requested interval
        spk0 = round((spks(rCh(ch), sw).tStamp(ii) - jj(1)) * 10);
        spk0((spk0 <= 0) | (spk0 > (size(tm, 1) * 10))) = [];
        spk1 = spks(rCh(ch), sw).ampl(ii);
        for t = 1:length(spk0)
            whtLFP(spk0(t), ch, z) = whtLFP(spk0(t), ch, z) + spk1(t);
        end
    end
end
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)
ajj = zeros(size(tm, 1), length(rCh));
for t = 1:size(whtLFP, 2)
    ajj(:, t) = resample(whtLFP(:, t), 100, 1000);%resampled lfp
end
whtLFP = ajj;

% % common reference
% wLFPcomm = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, 1:16, rawData, nlxVer);
% wLFPcomm = mean(wLFPcomm, 3); wLFPcomm = mean(wLFPcomm, 2);

%special channels
n = 0;%number of special channels
for t = 1:n
    whtLFP(:, t) = (whtLFP(:, t) - mean(whtLFP(:, t)));%amplification of signal (piezo)
    %whtLFP(:, t) = whtLFP(:, t) * max(max(abs(lfp(:, rCh((n + 1):end))))) / (0.3 * max(abs(lfp(:, rCh(t)))));
    whtLFP(:, t) = whtLFP(:, t) * mean(std(lfp(:, rCh((n + 1):end)), [], 1)) / (1 * std(lfp(:, rCh(t)), [], 1));
%     whtLFP(:, t) = ZavFilter(whtLFP(:, t), prm.Fs, 'high', 2, 2);%
end

%filtration
for ch = (n + 1):size(whtLFP, 2)
    %0
%     jj = 1:1e3:(2 * 60e3);
%     p1 = polyfit(tm(jj), whtLFP(jj, ch), 1);
%     whtLFP(:, ch) = whtLFP(:, ch) - (p1(1) * tm + p1(2));
    %1
%     whtLFP(:, ch) = locdetrend(whtLFP(:, ch), 1e3, [300, 10]);
    %2
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'bandpass', [1, 300], 2);%
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'low', 200, 2);%
%     whtLFP(:, ch) = locdetrend(whtLFP(:, ch), prm.Fs, [0.002, 0.001]);%for raw data (spikes)
    %3
%     whtLFP(:, ch) = ZavSpctSubtr(whtLFP(:, ch), []);%denoise by spectral subtraction
    
%     jj = ismember(-50:250, tm);
%     whtLFP(:, ch) = whtLFP(:, ch) - wLFPcomm(jj, rCh(ch));%common reference
%     ajj = mean(whtLFP(1:120, ch)); whtLFP(:, ch) = whtLFP(:, ch) - ajj;
end

if (diff(segmEdge) < 5e3)
    dCoeff(4) = 1; bufStr = 'ms';%time coefficient
elseif ((diff(segmEdge) >= 5e3) && (diff(segmEdge) < 5e4))
    dCoeff(4) = 1e3; bufStr = 's';%time coefficient
else
    dCoeff(4) = 60e3; bufStr = 'min';%time coefficient
end

clf, hold on
if pltCSD
    ajj = whtLFP - repmat(mean(whtLFP((tm >= 1) & (tm <= 3), :), 1), numel(tm), 1);%replace DC
%     ajj = ajj ./ repmat(mean(lfpVar(rCh, :), 2)' / min(mean(lfpVar(rCh, :), 2)), numel(tm), 1);%amplitude correction on impedance

%     ajj = filtfilt([0.607, 1, 0.607], 2.213, ajj')';%spatial filtration    
    lfpCSD = -diff(ajj, 2, 2);%CSD (current source density)
    k = 3;%round(diff(segmEdge) / 100) + 1;%interval of smoothing
    for t = 1:size(lfpCSD, 2) %run over (channels) depth
        lfpCSD(:, t) = smooth(lfpCSD(:, t), k);%smoothing along time
    end
    %lfpCSD = filtfilt([0.607, 1, 0.607], 2.213, lfpCSD')';%3-channel spatial smoothing as in elephy
    
    %interpolation (along channels)
    ajj = lfpCSD;
    spPrc = 0.1;%1;%spatial precision
    lfpCSD = zeros(size(ajj, 1), numel(1:spPrc:(numel(rCh) - 2)));
    for t = 1:size(lfpCSD, 1) %run over time
        lfpCSD(t, :) = interp1(1:(numel(rCh) - 2), ajj(t, :), 1:spPrc:(numel(rCh) - 2));%interpolation along channels
        %lfpCSD(t, :) = spline(1:(numel(rCh) - 2), ajj(t, :), 1:spprc:(numel(rCh) - 2));%spline interpolation along channels
    end
    
    imagesc(tm / dCoeff(4), 2:(numel(rCh) - 1), lfpCSD')
    set(gca, 'Color', 0.5 * ones(1, 3))%%gray rectangle background %imagesc(tm([1, end]) / dCoeff(4), [1, numel(rCh)], -caRng * ones(size(lfpCSD, 2), 2))
    colormap('default');
    hclrb = colorbar;%colobar handle
    set(get(hclrb, 'YLabel'), 'String', 'sink - source')

    caxis(caRng * [-1, 1])%colorbar range
    %caxis(max(max(abs(lfpCSD))) * [-1, 1]);%colorbar range
end

%%-plot LFP-%%
if (size(whtLFP, 1) > 1e4)
    ii = round(linspace(1, size(whtLFP, 1), 3e4));
else
    ii = 1:size(whtLFP, 1);
end
dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient
plot(tm(ii) / dCoeff(4), (whtLFP(ii, :) / dCoeff(3)) + repmat(1:numel(rCh), numel(ii), 1), 'k', 'LineWidth', 1)
%%=end of plot LFP=%%

% %%=plot MUA=%%
% if (numel(sn) == 1) %single sweep (segment)
%     jj = segms(sn, 1) + segmEdge;%absolute time of the segment boundary ( ms, e.i. from begin of file)
%     sw = segms(sn, 3);%sweep
%     for ch = 1:numel(rCh) %run over wanted channels
%         spkDens = spks(rCh(ch), sw).tStamp((spks(rCh(ch), sw).tStamp > jj(1)) & (spks(rCh(ch), sw).tStamp < jj(end))) - segms(sn, 1);%unites in the segment
%         plot([spkDens'; spkDens'] / dCoeff(4), ch + repmat(0.0 + [-0; 0.3], 1, numel(spkDens)), 'r', 'LineWidth', 1.5)
%     end
% end
% %%-end of plot MUA-%%

set(gca, 'YDir', 'reverse', 'XLim', segmEdge / dCoeff(4), 'YLim', [-0.05, numel(rCh) + 1.05])
set(get(gca, 'XLabel'), 'String', ['time, ', bufStr])

if isequal(hd.fFileSignature, 'eeg EDF+') %named labels
    bufStr = char;%2) label for man EEG
    for t = 1:numel(rCh)
        bufStr(t, 1:length(hd.label{rCh(t)})) = hd.label{rCh(t)};
    end
    z = segms(sn(end), 1) / (1e3 * 60);%absolute time (minutes)
    disp([num2str(floor(z)), 'min ', num2str(round((z - floor(z)) * 60)), 'sec'])%absolute time (minutes)
    if (numel(sn) == 1)
        [~, z] = min(abs(segms(sn, 1) - (hd.events.POS / zavp.rarStep)));
        disp([num2str(z), ' | ', hd.events.TYP{z}])
    end
else %common label for channels
    bufStr = num2str(rCh');%1) common label for channels
    disp([num2str((segms(sn([1, end]), 1)' + segmEdge) / 1e3), ' s'])
end

set(gca, 'YTick', 1:numel(rCh), 'YTickLabel', bufStr) %num2str(rCh'))%
set(get(gca, 'YLabel'), 'String', 'channels', 'Units', 'normalized', 'Position', [-0.06, 0.6, 1])
plot([0, 0], get(gca, 'YLim'), 'm', 'LineWidth', 2)%nul-time line

[~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
annotation('textarrow', 'Position', [0.05, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    'String', [num2str(round(dCoeff(2) * 10) / 10), ' \muV'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
else
    bufStr = [bufStr, num2str(sn)];
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if (isempty(p2) || (p2 < p1)), p2 = length(zavp.file) - 1; end
if pltCSD
    fN = [zavp.file(p1:p2), bufStr, '_CSD'];
else
    fN = [zavp.file(p1:p2), bufStr, '_lfp'];
end
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title(fN)


%% === multi-dimensional LFP, CSD, MUA === %%
rawData = 0;%false;%if raw data needed then rawData = 1(true)
pltCSD = 1;%plot CSD?
dCoeff(1) = -1.0; dCoeff(2) = 200;%max(max(abs(whtLFP)));%max(max(abs(lfp(:, rCh((n + 1):end)))));%characteristic sizes (amplitude coefficient)
caRng = 100;%max(max(abs(lfpCSD)));%range of colormap for CSD

segmEdge = [-250, 250];%left and right shifts from synchro-point (ms)
rCh = 1:32;%channels to be read and draw (1:hd.nADCNumChannels)

% %poly3
% probX = [[1:10; 1 * ones(1, 10)], [11:22; 1.5 * ones(1, 12)], [23:32; 2 * ones(1, 10)]];%x-coordinates of probe sites
% probY = [[1:10; 2.5:11.5], [11:22; 1:12], [23:32; 2.5:11.5]];%y-coordinates of probe sites

%A4x8 (200x50)
probX = [[1:8; 1 * ones(1, 8)], [9:16; 5 * ones(1, 8)], [17:24; 9 * ones(1, 8)], [25:32; 13 * ones(1, 8)]];
probY = [[1:8; 1:8], [9:16; 1:8], [17:24; 1:8], [25:32; 1:8]];

hx = min(diff(unique(probX(2, :))));%spacing over x-coordinate
hy = min(diff(unique(probY(2, :))));%spacing over y-coordinate
[xgrd, ygrd] = meshgrid(min(probX(2, :)):hx:max(probX(2, :)), min(probY(2, :)):hy:max(probY(2, :)));

sn = 1:size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));
if isempty(sn), return; end

% whtLFP = ZavFilter(lfp, prm.Fs, 'high', 2, 2); whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, whtLFP, rCh, rawData, nlxVer);%lfp phased with respect to synchro-events
whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to synchro-events
whtLFP = whtLFP - repmat(whtLFP(max(5, 1), :, :), size(whtLFP, 1), 1);%DC-shift
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)

% %common reference
% wLFPcomm = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, 1:16, rawData, nlxVer);
% wLFPcomm = mean(wLFPcomm, 3); wLFPcomm = mean(wLFPcomm, 2);


if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k))' / k;%time matrix for raw or downsampled data (ms)

%filtration
for ch = 1:size(whtLFP, 2)
    %1
%     whtLFP(:, ch) = locdetrend(whtLFP(:, ch), 1e3, [100, 5]);
    %2
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'high', 8, 2);%
%     whtLFP(:, ch) = ZavFilter(whtLFP(:, ch), prm.Fs, 'low', 30, 2);%
    %3
%     whtLFP(:, ch) = ZavSpctSubtr(whtLFP(:, ch), []);%denoise by spectral subtraction
    
%     whtLFP(:, ch) = whtLFP(:, ch) - wLFPcomm;%common reference
%     ajj = mean(whtLFP(:, ch)); whtLFP(:, ch) = whtLFP(:, ch) - ajj;
end

clf, hold on
if (diff(segmEdge) > 5e3), dCoeff(4) = 1e3; bufStr = 's'; else dCoeff(4) = 1; bufStr = 'ms'; end %time coefficient
if pltCSD
    ajj = whtLFP - repmat(mean(whtLFP(1:20, :), 1), numel(tm), 1);%replace DC
    %ajj = locdetrend(whtLFP, 1e3, [0.3, 0.1]);
    ajj = ajj ./ repmat(mean(lfpVar(rCh, :), 2)' / min(mean(lfpVar(rCh, :), 2)), numel(tm), 1);%amplitude correction

%     %lfpCSD = -diff(ZavFilter(ajj, prm.Fs, 'low', 150, 2), 2, 2);%CSD (current source density)
%     for t = 1:size(ajj, 2) %run over time
%         ajjG = griddata(probX, probY, ajj(:, t), xgrd, ygrd);
%     end
%     lfpCSD = -del2(ajjG, hx, hy);%CSD (current source density)

    ajjG = griddata(probX(2, :), probY(2, :), whtLFP(t, 1:32), xgrd, ygrd);
    lfpCSD = -del2(ajjG, hx, hy);
    contour(xgrd, ygrd, lfpCSD), set(gca, 'YDir', 'reverse'), caxis([-0.1, 0.1])

    k = 3;%round(diff(segmEdge) / 100) + 1;%interval of smoothing
    for t = 1:size(lfpCSD, 2)
        lfpCSD(:, t) = smooth(lfpCSD(:, t), k);%smoothing
    end
    %interpolation (along channels)
    ajj = lfpCSD;
    spPrc = 0.1;%1;%spatial precision
    lfpCSD = zeros(size(ajj, 1), numel(1:spPrc:(numel(rCh) - 2)));
    for t = 1:size(lfpCSD, 1)
        lfpCSD(t, :) = interp1(1:(numel(rCh) - 2), ajj(t, :), 1:spPrc:(numel(rCh) - 2));%interpolation along channels
    end

    imagesc(tm / dCoeff(4), 2:(numel(rCh) - 1), lfpCSD'), set(gca, 'YDir', 'reverse')%normal
    set(gca, 'Color', 0.5 * ones(1, 3))
    colormap('default');
    hclrb = colorbar;%colobar handle
    set(get(hclrb, 'YLabel'), 'String', 'sink - source')

    caxis(caRng * [-1, 1])%colorbar range
end

%plot spikes and lfp
if (numel(sn) == 1) %single sweep (segment)
    jj = segms(sn, 1) + segmEdge;%absolute time of the segment boundary ( ms, e.i. from begin of file)
    sw = segms(sn, 3);%sweep
    for ch = 1:numel(rCh) %run over wanted channels
        spkDens = spks(rCh(ch), sw).tStamp((spks(rCh(ch), sw).tStamp > jj(1)) & (spks(rCh(ch), sw).tStamp < jj(end))) - segms(sn, 1);%unites in the segment
        plot([spkDens'; spkDens'] / dCoeff(4), ch - repmat([0.1; 0.4], 1, numel(spkDens)), 'r', 'LineWidth', 1.5)
        %plot([spkDens'; spkDens'] / dCoeff(4), ch - repmat([0.1; 0.4], 1, numel(spkDens)) + 0.03 * repmat(rand(1, numel(spkDens)), 2, 1), 'r', 'LineWidth', 1.5)
    end
end
dCoeff(3) = dCoeff(1) * dCoeff(2);%scaling coefficient
plot(tm / dCoeff(4), (whtLFP / dCoeff(3)) + repmat(1:numel(rCh), numel(tm), 1), 'k', 'LineWidth', 1.5)

set(gca, 'YDir', 'reverse', 'XLim', segmEdge / dCoeff(4), 'YLim', [-0.05, numel(rCh) + 1.05])
set(get(gca, 'XLabel'), 'String', ['time, ', bufStr])

if isequal(hd.fFileSignature, 'eeg EDF+') %named labels
    bufStr = char;%2) label for man EEG
    for t = 1:numel(rCh)
        bufStr(t, 1:length(hd.label{rCh(t)})) = hd.label{rCh(t)};
    end
    z = segms(sn(end), 1) / (1e3 * 60);%absolute time (minutes)
    disp([num2str(floor(z)), 'min ', num2str(round((z - floor(z)) * 60)), 'sec'])%absolute time (minutes)
    if (numel(sn) == 1)
        [~, z] = min(abs(segms(sn, 1) - (hd.events.POS / zavp.rarStep)));
        disp([num2str(z), ' | ', hd.events.TYP{z}])
    end
else %common label for channels
    bufStr = num2str(rCh');%1) common label for channels
    disp([num2str((segms(sn([1, end]), 1)' + segmEdge) / 1e3), ' s'])
end

set(gca, 'YTick', 1:numel(rCh), 'YTickLabel', bufStr) %num2str(rCh'))%
set(get(gca, 'YLabel'), 'String', 'channels', 'Units', 'normalized', 'Position', [-0.06, 0.6, 1])
plot([0, 0], get(gca, 'YLim'), 'r', 'LineWidth', 2)%nul-time line

[~, p1] = ds2nfu([0, 0], [0, abs(1 / dCoeff(1))]);%diff(p1) = length of scale bar in normalized units
annotation('textarrow', 'Position', [0.05, 0.15 0, diff(p1)], 'LineWidth', 5, 'Units', 'normalized', 'TextRotation', 90, 'HeadStyle', 'none',...
    'String', [num2str(round(dCoeff(2) * 10) / 10), ' \muV'], 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'left');

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
else
    bufStr = [bufStr, num2str(sn)];
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if (isempty(p2) || (p2 < p1)), p2 = length(zavp.file) - 1; end
if pltCSD
    fN = [zavp.file(p1:p2), bufStr, '_CSD'];
else
    fN = [zavp.file(p1:p2), bufStr, '_lfp'];
end
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title(fN)

%% === power spectrum === %%
rawData = 0;%if raw data needed then rawData = 1(true)

segmEdge = [-4000, 0];%left and right shifts from synchro-point (ms)
rCh = 3;%channels to be read and draw
sn = 580:680;%690 + (0:50);%:size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
whtLFP = whtLFP - repmat(median(whtLFP(1:min(10, size(whtLFP, 1)), :, :), 1), size(whtLFP, 1), 1);%DC-shift
whtLFP = squeeze(whtLFP);

if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
prm.fpass = [1, 100];
prm.tapers = [3, 5];%tapers %[diff(prm.fpass), diff(segmEdge) * 1e-3, 1]
prm.pad = 2;%padding factors
prm.trialave = 1;%average over trials/channels when 1

% %get new segments from wavelets
% segms1 = segms(sn, :);
% [wavlCoef, wavlFrq] = Wavelet(whtLFP, prm.fpass, prm.Fs, 6, 0);%Scalogram, scalogram(from SpecM)
% wavlCoef = abs(wavlCoef) .^ 2;%mean(abs(wavlCoef) .^ 2, 3);
% for t = 1:numel(sn)
%     pwrSpct = sum(wavlCoef(:, :, t), 2);
%     [~, z] = max(pwrSpct);
%     segms1(t, 1) = segms1(t, 1) + z + segmEdge(1);
% end
% segmEdge = [-750, 750];%left and right shifts from synchro-point (ms)
% whtLFP = ZavSynchLFP(zavp, hd, segms1, segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
% whtLFP = squeeze(whtLFP);


% %1
% arOrder = 1;%autoregression model order
% arModel = zeros(length(sw), arOrder + 2);%coefficients of autoregression model
% % whtLFP = ZavFilter(whtLFP, prm.Fs, 'low', 150, 2);
% for t = 1:numel(sn)
%     [~, arModel(t, :)] = WhitenSignal(whtLFP(:, t), [], [], [], arOrder);%WhitenSignal(squeeze(lfpShft(1:rspLen, ch, t)), [], [], [], arOrder);%autoregressive mode
% end
% arModel = mean(arModel, 1);
% whtLFP = WhitenSignal(whtLFP, [], [], arModel, arOrder);%
% 
% %2
% whtLFP = ZavFilter(whtLFP, prm.Fs, 'high', 200, 2);
% whtLFP = locdetrend(whtLFP, prm.Fs, [0.1, 0.01]);%

% %3
% for t = 1:size(whtLFP, 2)
%     whtLFP(:, t) = whtLFP(:, t) - max(whtLFP(:, t));
% %     whtLFP = ZavFilter(whtLFP, prm.Fs, 'high', 2, 2);
%     whtLFP(:, t) = ZavSpctSubtr(whtLFP(:, t), []);%denoise by spectral subtraction
% end

[pwrSpct, frq, err] = mtspectrumc(whtLFP, prm);%mean power spectrum (mean by sweeps)
clf
% plot(log(frq), log(pwrSpct), 'LineWidth', 3)%double-log axis
% [~, jj] = ZavPowerAppr(frq', pwrSpct);%spectrum subtraction
% pwrSpct = pwrSpct - jj;
err(1, :) = smooth(err(1, :), 5); err(2, :) = smooth(err(2, :), 5);
plot(repmat(frq, 2, 1), [err(1, :); err(2, :)], 'r', 'LineWidth', 2); hold on
plot(frq, smooth(pwrSpct, 3), 'k', 'LineWidth', 2)

set(gca, 'XLim', [0.85 * prm.fpass(1), 1.02 * prm.fpass(2)]);
xlim(frq([1, end])), %ylim([0 1]) %set(gca, 'XTick', (0:20:prm.fpass(2)))
set(get(gca, 'XLabel'), 'String', 'frequency, Hz')
set(get(gca, 'YLabel'), 'String', 'power, \muV^2/Hz')
set(gca, 'YScale', 'log')
% set(gca, 'XScale', 'log')

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
else
    bufStr = [bufStr, num2str(sn)];
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if isempty(p2), p2 = length(zavp.file) - 1; end
fN = [zavp.file(p1:p2), bufStr, '_pwrSpct'];
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title(fN)


%% === spectrogam === %%
%fLFP = ZavFilter(lfp(:, 17:32), 1e3, 'high', 1, 2);
rawData = false;%if raw data needed then rawData = 1(true)
segmEdge = [1000, 2000];%left and right shifts from synchro-point (ms)
rCh = 13;%[2, 7:2:13];%channels to be read and draw
sn = 178:323;%1:size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));
if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
prm.fpass = [2, 80];
prm.tapers = [3, 5];%tapers
prm.pad = 2;%padding factors
prm.trialave = 0;%average over trials/channels when 1
[~, frq] = mtspectrumc(randn(diff(segmEdge) + 1, 1), zavp.prm);%power spectrum
m = 3;
pwrSpct = Inf(length(frq), numel(sn) * m);
for zT = 1:m
    segmEdge = segmEdge + (zT - 1) * 1000;
    whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, 1:16, rawData, nlxVer);
    wLFPcomm = ZavFilter(squeeze(mean(whtLFP, 2)), prm.Fs, 'bandpass', [47, 53], 1);
    
    whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, fLFP, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
    whtLFP = squeeze(whtLFP) - wLFPcomm;
    
    %whtLFP = ZavFilter(whtLFP, prm.Fs, 'stop', [46, 57], 1);%
    for t = 1:numel(sn) %run over segmens
        pwrSpct(:, zT + ((t - 1) * m)) = mtspectrumc(whtLFP(:, t), zavp.prm);%power spectrum
    end
end
for t = 1:size(pwrSpct, 1)
    pwrSpct(t, :) = smooth(pwrSpct(t, :), 3);
end
for t = 1:size(pwrSpct, 2)
    pwrSpct(:, t) = smooth(pwrSpct(:, t), 3);
end
tm = sort([segms(sn, 1); (segms(sn, 1) + 1000); (segms(sn, 1) + 2000)])  - segms(298, 1);
if (diff(tm([1, end])) < 5e3)
    dCoeff(4) = 1; bufStr = 'ms';%time coefficient
elseif ((diff(tm([1, end])) >= 5e3) && (diff(segmEdge) < 5e4))
    dCoeff(4) = 1e3; bufStr = 's';%time coefficient
else
    dCoeff(4) = 60e3; bufStr = 'min';%time coefficient
end
clf
imagesc(tm / dCoeff(4), frq, log(pwrSpct)), set(gca, 'YDir', 'normal'), hold on
colormap('default'); colorbar;
plot([0, 0], get(gca, 'YLim'), 'm', 'LineWidth', 2)%stim-line
set(get(gca, 'XLabel'), 'String', ['time, ', bufStr])
set(get(gca, 'YLabel'), 'String', 'Frequency, Hz')
ylim(prm.fpass)%([5, 30])
caxis([0, 5])
set(gca, 'YTick', 10:20:80)


%% === spectral deep profile === %%
rawData = 0;%if raw data needed then rawData = 1(true)

segmEdge = [50, 200];%left and right shifts from synchro-point (ms)
rCh = 17:32;%channels to be read and draw
sn = 500:size(segms, 1);%wanted sweeps or segments
sn = sn(isfinite(segms(sn, 1)));

if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
prm.fpass = [8, 100];
prm.tapers = [3, 5];%tapers %[diff(prm.fpass), diff(segmEdge) * 1e-3, 1]
prm.pad = 2;%padding factors
prm.trialave = 1;%average over trials/channels when 1

    
[pwrSpct, frq] = mtspectrumc(rand(diff(segmEdge) + 1, 1), zavp.prm);
pwrDeep = zeros(length(pwrSpct), numel(rCh));
for ch = 1:numel(rCh)
    whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh(ch), rawData, nlxVer);%lfp phased with respect to stimulus moments
    whtLFP = whtLFP - repmat(median(whtLFP(1:min(10, size(whtLFP, 1)), :, :), 1), size(whtLFP, 1), 1);%DC-shift
    whtLFP = squeeze(whtLFP);

%     %1
%     arOrder = 1;%autoregression model order
%     arModel = zeros(length(sw), arOrder + 2);%coefficients of autoregression model
%     % whtLFP = ZavFilter(whtLFP, prm.Fs, 'low', 150, 2);
%     for t = 1:numel(sn)
%         [~, arModel(t, :)] = WhitenSignal(whtLFP(:, t), [], [], [], arOrder);%WhitenSignal(squeeze(lfpShft(1:rspLen, ch, t)), [], [], [], arOrder);%autoregressive mode
%     end
%     arModel = mean(arModel, 1);
%     whtLFP = WhitenSignal(whtLFP, [], [], arModel, arOrder);%
% 
%     %2
%     whtLFP = ZavFilter(whtLFP, prm.Fs, 'high', 200, 2);
%     whtLFP = locdetrend(whtLFP, prm.Fs, [0.1, 0.01]);%
% 
%     %3
%     for t = 1:size(whtLFP, 2)
%         whtLFP(:, t) = whtLFP(:, t) - max(whtLFP(:, t));
%     %     whtLFP = ZavFilter(whtLFP, prm.Fs, 'high', 2, 2);
%         whtLFP(:, t) = ZavSpctSubtr(whtLFP(:, t), []);%denoise by spectral subtraction
%     end

    [pwrSpct, frq] = mtspectrumc(whtLFP, zavp.prm);%mean power spectrum (mean by sweeps)
    pwrDeep(:, ch) = pwrSpct;
end

clf
imagesc(frq, rCh, log(pwrDeep)'), colorbar, colormap('default')

xlim(frq([1, end])), ylim(rCh([1, end])) %set(gca, 'XTick', (0:20:prm.fpass(2)))
set(get(gca, 'XLabel'), 'String', 'frequency, Hz')
set(get(gca, 'YLabel'), 'String', 'channels')

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
else
    bufStr = [bufStr, num2str(sn)];
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if isempty(p2), p2 = length(zavp.file) - 1; end
fN = [zavp.file(p1:p2), bufStr, '_deep-spectra'];
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title(fN)



%% === wavelet transformation === %%
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [0, 240e3] + 0e3;%left and right shifts from synchro-point (ms)
rCh = 5;%channels to be read and draw
sn = 51;%:size(segms, 1);%wanted sweeps or segments

for ch = 1:numel(rCh)
    subplot(numel(rCh), 1, ch)
whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh(ch), rawData, nlxVer);%lfp phased with respect to stimulus moments
whtLFP = squeeze(whtLFP);%average by sweeps (segments)
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)

prm.fpass = [6, 16];
if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data

% %1
% arOrder = 1;%autoregression model order
% arModel = zeros(length(sw), arOrder + 2);%coefficients of autoregression model
% whtLFP = ZavFilter(whtLFP, prm.Fs, 'low', 150, 2);
% for t = 1:numel(sn)
%     [~, arModel(t, :)] = WhitenSignal(whtLFP(:, t), [], [], [], arOrder);%WhitenSignal(squeeze(lfpShft(1:rspLen, rCh, t)), [], [], [], arOrder);%autoregressive mode
% end
% arModel = mean(arModel, 1);
% whtLFP = WhitenSignal(whtLFP, [], [], arModel, arOrder);%

%2
% whtLFP = ZavFilter(whtLFP, prm.Fs, 'high', 1, 1);

% %3
% for t = 1:size(whtLFP, 2)
%     whtLFP(:, t) = whtLFP(:, t) - max(whtLFP(:, t));
%     whtLFP = ZavFilter(whtLFP, prm.Fs, 'high', prm.fpass(1), 2);
%     whtLFP(:, t) = ZavSpctSubtr(whtLFP(:, t), []);%denoise by spectral subtraction
% end

[wavlCoef, wavlFrq] = Wavelet(whtLFP, prm.fpass, prm.Fs, 6, 0);%Scalogram, scalogram(from SpecM)
wavlCoef = abs(wavlCoef .^ 2);%mean(abs(wavlCoef) .^ 2, 3);
wavlCoef = mean(wavlCoef, 3);%mean by segments
tm = tm(1:5:end);
wavlCoef = wavlCoef(1:5:end, :);
for t = 1:size(wavlCoef, 2) %run over frequency
    wavlCoef(:, t) = smooth(wavlCoef(:, t), 255);
end
% for t = 1:size(wavlCoef, 1)
%     wavlCoef(t, :) = log(wavlCoef(t, :));
%     wavlCoef(t, wavlCoef(t, :) < 0) = 0;
% end
wavlCoef(wavlCoef < 4.6e4) = wavlCoef(wavlCoef < 4.6e4) / 1.2;

clf
if (diff(segmEdge) > 5e3), dCoeff(4) = 1e3; bufStr = 's'; else dCoeff(4) = 1; bufStr = 'ms'; end %time coefficient
imagesc(tm / dCoeff(4), wavlFrq, wavlCoef'), set(gca, 'YDir', 'normal'), hold on
colormap('default');
% colorbar;
% clr_mp = colormap; clr_mp(1, :) = 1; colormap(clr_mp);
plot([0, 0], get(gca, 'YLim'), 'm', 'LineWidth', 2)%stim-line
set(get(gca, 'XLabel'), 'String', ['time, ', bufStr])
set(get(gca, 'YLabel'), 'String', 'Frequency, Hz')

ylim(prm.fpass)%([5, 30])
caxis([0, 7e4])%color range

% wvName = 'morl';%name of wavelet
% f_centr = m_centfrq(wvName);%frequency of central wavelet splash (peak)
% frq = prm.fpass(1):0.5:prm.fpass(2);%frequency scale
% wavlScl = (f_centr * prm.Fs) ./ frq;%wavelet scale
% wavlFrq = m_scal2frq(wavlScl, wvName, 1 / prm.Fs);%pseudo-frequencies            
% wavlCoef = zeros(length(wavlScl), size(whtLFP, 1));%
% 
% for t = 1:size(whtLFP, 2)
%     wavlCoef = wavlCoef + abs(m_cwt(whtLFP(:, t), wavlScl, wvName));
% end
% k = round(size(wavlCoef, 2) / 200);
% for t = 1:size(wavlScl, 2)
%     wavlCoef(:, t) = smooth(wavlCoef(:, t), k);%smoothing wavelet coefficients
% end
% for t = 1:size(wavlScl, 1)
%     wavlCoef(t, :) = smooth(wavlCoef(t, :), 50);%smoothing wavelet coefficients
% end
% wavlCoef = wavlCoef / size(whtLFP, 2);
% %wavlCoef(:, end + 1) = 0;%additional line for visibility
% 
% wavlCoef(end, [1, end]) = [1e-10, 1e6];%range scale 12668
% wavlCoef = log(wavlCoef);%log scale
% %wavlFrq = log(wavlFrq);%log-log-scale
% %wavlCoef(end, end) = 8;%colorbar range
% %80;%whiten %8;%log %2180;%linear
% 
% clf
% imagesc(tm, wavlFrq, wavlCoef), set(gca, 'YDir', 'normal')
% colormap('default'); colorbar; hold on
% plot([0, 0], get(gca, 'YLim'), 'r', 'LineWidth', 2)%stim-line
% set(gca, 'XLim', segmEdge)
% set(get(gca, 'XLabel'), 'String', 'time, ms')
% set(get(gca, 'YLabel'), 'String', 'frequency, Hz')
% 
% %nonlinear colormap
% colormap('default'); clrMap = colormap;
% t = 1:size(clrMap, 1);
% indx = round(64 * ((atan((t - 50) / 1.3) + (pi / 2)) / pi));%sigmoid
% clrMap = clrMap(indx, :); colormap(clrMap)
% 
% % [wCf, f, tm] = Wavelet(whtLFP, [1, 120], prm.Fs, 6, 1);        
% % wavlCoef = mean(abs(wCf) .^ 2, 3);
% % pcolor(tm, f, log(wavlCoef)'), shading flat;
% % % Scalogram(whtLFP(:, 3));
end

if isequal(hd.fFileSignature, 'eeg EDF+') %named labels
    z = segms(sn(end), 1) / (1e3 * 60);%absolute time (minutes)
    disp([num2str(floor(z)), 'min ', num2str(round((z - floor(z)) * 60)), 'sec'])%absolute time (minutes)
else %common label for channels
    bufStr = num2str(rCh');%1) common label for channels
    disp([num2str((segms(sn([1, end]), 1)' + segmEdge) / 1e3), ' s'])
end

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
else
    bufStr = [bufStr, num2str(sn)];
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if (isempty(p2) || (p2 < p1)), p2 = length(zavp.file) - 1; end
fN = [zavp.file(p1:p2), bufStr, '_wavLet'];
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
title(fN)


%% === MUA-LFP coherence === %%
rawData = 0;%if raw data needed then rawData = 1(true)

segmEdge = [0, 500];%left and right shifts from synchro-point (ms)
rCh = 27;%channels to be read and draw
sn = 70:78;%1:size(segms, 1);%wanted sweeps or segments

sn = sn(isfinite(segms(sn, 1)));
if isempty(sn), return; end

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to synchro-events
whtLFP = whtLFP - repmat(median(whtLFP(1:min(10, size(whtLFP, 1)), :, :), 1), size(whtLFP, 1), 1);%DC-shift
whtLFP = squeeze(whtLFP);%

if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
prm.fpass = [15, 120];
prm.tapers = [3, 5];%tapers
prm.pad = 2;%padding factor
prm.trialave = 1;%average

% %1
% arOrder = 1;%autoregression model order
% arModel = zeros(length(sw), arOrder + 2);%coefficients of autoregression model
% % whtLFP = ZavFilter(whtLFP, prm.Fs, 'low', 150, 2);
% for t = 1:numel(sn)
%     [~, arModel(t, :)] = WhitenSignal(whtLFP(:, t), [], [], [], arOrder);%WhitenSignal(squeeze(lfpShft(1:rspLen, ch, t)), [], [], [], arOrder);%autoregressive mode
% end
% arModel = mean(arModel, 1);
% whtLFP = WhitenSignal(whtLFP, [], [], arModel, arOrder);%

% %2
% whtLFP = ZavFilter(whtLFP, prm.Fs, 'bandpass', [200, 800], 2);% mode 4 - Butterworth filter
% whtLFP = locdetrend(whtLFP, prm.Fs, [0.1, 0.01]);%

% %3
% for t = 1:size(whtLFP, 2)
%     whtLFP(:, t) = whtLFP(:, t) - max(whtLFP(:, t));
%     whtLFP(:, t) = ZavSpctSubtr(whtLFP(:, t), []);%denoise by spectral subtraction
% end


clear spkNew spkBef
spkNew(1:numel(sn)) = struct('s', []);%spikes in seconds
spkBef(1:numel(sn)) = struct('s', []);%spikes in seconds

% allSpk = [];%all spikes
% tW = segms(t, 1);%cutting window
% for t = sn
%     allSpk = [allSpk; spks(ch, t).tStamp(spks(ch,t).tStamp > tW)];
% end
% [spkDens, bins] = hist(allSpk, 100); spkDens = smooth(spkDens, 5);
% k = find(spkDens >= (0.2 * max(spkDens)), 1, 'first');%

z = 1;
for t = sn
    tW = segms(t, 1) + segmEdge;%cutting window (signal after beginning)
    sw = segms(t, 3);%number of sweep
    spkNew(z).s = spks(rCh(1), sw).tStamp((spks(rCh(1), sw).tStamp >= tW(1)) & (spks(rCh(1), sw).tStamp <= tW(2))) / 1e3;%spikes time (seconds)
    if ~isempty(spkNew(z).s)
        spkNew(z).s = spkNew(z).s - (tW(1) / 1e3);%putting to zero
    end
    tW = segms(t, 1) - segmEdge(end:-1:1);%signal befor beginning
    spkBef(z).s = spks(rCh(1), sw).tStamp((spks(rCh(1), sw).tStamp >= tW(1)) & (spks(rCh(1), sw).tStamp <= tW(2))) / 1e3;%spikes time (seconds)
    if ~isempty(spkBef(z).s)
        spkBef(z).s = spkBef(z).s - spkBef(z).s(1);%putting to zero
    end
    z = z + 1;
end
%[chr, phi, S12, S1, S2, frq, zerosp, confC, phistd, sCerr] = coherencycpt() %chronux %%jackknife(x)
[chr, ~, ~, ~, ~, frq, ~, ~, ~, err] = coherencycpt(whtLFP, spkNew, zavp.prm);
%[chr, ~, ~, S1, S2, ~, ~, ~, ~, err, J1, J2] = ZavCoherencycpt(whtLFP, spkNew, zavp.prm);
%err = [chr - err(1, :)', err(2, :)' - chr];%errors (not errorbars as in coherencycpt)

% whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), -segmEdge(end:-1:1), lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
% whtLFP = squeeze(whtLFP);
% [chrB, ~, ~, ~, ~, frq, ~, ~, ~, sCerrB, J1b, J2b] = coherencycpt(whtLFP, spkBef, zavp.prm);%signal befor stimulus
% %[chrb, ~, ~, ~, ~, frq, ~, ~, ~, sCerrB, J1b, J2b] = ZavCoherencycpt(whtLFP, spkBef, zavp.prm);%signal befor stimulus
% %sCerrB = [chrB - sCerrB(1, :)', sCerrB(2, :)' - chrB];%errors (not errorbars as in coherencycpt)

clf
% subplot(2, 1, 1)
plot(repmat(frq, 2, 1), [err(1, :); err(2, :)], 'r', 'LineWidth', 2); hold on
plot(frq, chr, 'k', 'LineWidth', 1.5)
% plot(repmat(frq, 2, 1), [sCerrB(1, :); sCerrB(2, :)], 'g', 'LineWidth', 2)
% plot(frq, chrB, 'g', 'LineWidth', 1.5)%signal befor stimulus
% set(gca, 'YLim', [0, max([err(1, :), 1])])%[min(sCerrB(1, :)), max([err(1, :), 1])])

p1 = get(gca, 'YLim'); p2 = get(gca, 'XLim');
set(gca, 'YLim', p1 + diff(p1) * [-0.03, 0.03], 'XLim', p2 + diff(p2) * [-0.03, 0.03])
set(gca, 'XTick', (0:20:prm.fpass(2)))
set(get(gca, 'YLabel'), 'String', 'coherence')

% J1 = reshape(J1, size(J1, 1), size(J1, 2) * size(J1, 3));%J1 = squeeze(mean(J1, 2));%
% J2 = reshape(J2, size(J2, 1), size(J2, 2) * size(J2, 3));%J2 = squeeze(mean(J2, 2));%
% J1b = reshape(J1b, size(J1b, 1), size(J1b, 2) * size(J1b, 3));%J1b = squeeze(mean(J1b, 2));%
% J2b = reshape(J2b, size(J2b, 1), size(J2b, 2) * size(J2b, 3));%J2b = squeeze(mean(J2b, 2));%
% [dz, vdz, Adz] = two_group_test_coherence(J1, J2, J1b, J2b, 0.01);
% t = find(Adz < 1, 1, 'last');
% if ~isempty(t)
%     plot(frq(t) * [1, 1], 1.3 * p1, 'k', 'LineWidth', 3)
% end

% subplot(2, 1, 2)
% errorbar(frq, dz, vdz, vdz, 'b', 'LineWidth', 1.5), hold on, set(gca, 'YLim', [min(dz - vdz),  max(dz + vdz)])%[0, 20])%
% p1 = get(gca, 'YLim'); %p2 = get(gca, 'XLim');
% set(gca, 'YLim', p1 + diff(p1) * [-0.06, 0.03], 'XLim', p2 + diff(p2) * [-0.03, 0.03])
% set(gca, 'XTick', (0:20:prm.fpass(2)))
% set(get(gca, 'YLabel'), 'String', 'dz')
% set(get(gca, 'XLabel'), 'String', 'frequency, Hz')
% if ~isempty(t)
%     plot(frq(t) * [1, 1], 1.3 * p1, 'k', 'LineWidth', 3)
% end

ylim([0, 1])
xlim(prm.fpass)

bufStr = '_s';%litter for define number of sweep or segment
if (numel(sn) > 1)
    bufStr = [bufStr, num2str(sn(1)), '-', num2str(sn(end))];
else
    bufStr = [bufStr, num2str(sn)];
end
p1 = find(zavp.file(1:end - 1) == '\', 1, 'last') + 1;%
p2 = find(zavp.file == '.', 1, 'last') - 1; if isempty(p2), p2 = length(zavp.file) - 1; end
fN = [zavp.file(p1:p2), bufStr, '_mua-lfp_cohr'];
set(gcf, 'FileName', ['c:\Users\Zav\Pictures\', fN, '.emf'], 'Name', fN, 'NumberTitle', 'off')
fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
% subplot(2, 1, 1)
title(fN)


%% === clusterization: multichannel MUA === %%
rCh = 2:15;
lL = round(0.7 / (zavp.siS * 1e3));%left lenght (1ms)
rL = round(1 / (zavp.siS * 1e3));%right lenght (1ms)

% %loading data
% clear cmua
% cmua(1:max(rCh)) = struct('ch', Inf);%, 's', []);
% nlxVer = 0;%old format of Cheetah log-file
% data = ZavLoadData(zavp.file, hd, (rCh(1) - 1):rCh(1), 'a', 0, 'e', nlxVer);%read one channel all sweeps
% data = ZavFilter(data, 1 / zavp.siS, 'high', 300, 2);%highpass filtration
% data = ZavFilter(data, 1 / zavp.siS, 'low', 5000, 2);%lowpass filtration
% for ch = rCh %run over non-edge LFP-channels
%     cmua(ch).ch = ch;
%     data(:, 3) = ZavLoadData(zavp.file, hd, ch + 1, 'a', 0, 'e', nlxVer);%read one channel all sweeps
%     data(:, 3) = ZavFilter(data(:, 3), 1 / zavp.siS, 'high', 300, 2);%highpass filtration
%     data(:, 3) = ZavFilter(data(:, 3), 1 / zavp.siS, 'low', 5000, 2);%lowpass filtration
%     
%     jj = spks(ch, 1).ampl < (-6 * lfpVar(ch));%large spikes
%     spkDens = (spks(ch, 1).tStamp(jj) * zavp.rarStep) + 1;%moments of spikes (samples, from 1)
%     whtLFP = spks(ch, 1).ampl(jj);
%     ii = (spkDens < (2 * lL)) | (spkDens > (size(data, 1) - (2 * rL)));
%     spkDens(ii) = [];
%     whtLFP(ii) = [];
%     spkDens = round(spkDens);
%     
%     %delete closest neighbors
%     m = 0.6 * zavp.rarStep;%0.6ms - minimal distance between spikes
%     t = 1;
%     while (t < length(spkDens))
%         if ((spkDens(t + 1) - spkDens(t)) < m)
%             if (whtLFP(t + 1) > whtLFP(t))
%                 spkDens(t + 1) = [];
%                 whtLFP(t + 1) = [];
%             else
%                 spkDens(t) = [];
%                 whtLFP(t) = [];
%             end
%         else
%             t = t + 1;
%         end
%     end
%     spkN = numel(spkDens);%spikes number
%     if (spkN > 0)
%         cmua(ch).s(1:spkN) = struct('tS', Inf, 'prm', Inf(1, 6), 'wave', zeros(lL + rL + 1, 3), 'clst', Inf);
%         for z = 1:spkN %run over unites
%             cmua(ch).s(z).tS = spkDens(z);%moment of current spike (samples, from 1)
%             cmua(ch).s(z).prm = min(data(spkDens(z) + (-3:3), [2, 1, 3]), [], 1);%amplitudes of all spike-channel (3 channels)
% 
%             p1 = ZavFindMins(-data(spkDens(z) + (-lL:0), 2));
%             if isempty(p1), p1 = 1; end
%             cmua(ch).s(z).prm(4) = data(spkDens(z) + p1(end) - lL - 1, 2);%left maximum of central channel
%             p2 = ZavFindMins(-data(spkDens(z) + (0:rL), 2));
%             if isempty(p2), p2 = rL + 1; end
%             cmua(ch).s(z).prm(5) = data(spkDens(z) + p2(1) - 1, 2);%right maximum of central channel
%             cmua(ch).s(z).prm(6) = p2(1) + (lL - p1(end) + 1);%duration
%             cmua(ch).s(z).wave = data(spkDens(z) + (-lL:rL), :);%waveforms
%         end
%     else
%         cmua(ch).s = [];
%     end
%     data(:, 1) = data(:, 2);
%     data(:, 2) = data(:, 3);
% end
% %save clstXXXX-XX-XX_XX-XX-XX cmua
% %end of loading data

%clustering
clc
%matlabpool('open', 5);%call to open the distributed processing
for n = 1:length(cmua) %run over non-edge LFP-channels
    spkN = length(cmua(n).s);%spikes number
    if (spkN > 0) %spikes was detected
        %1)PCA
        err = zeros(spkN, (lL + rL + 1) * 3);
        for t = 1:spkN %run over spikes of current channel
            err(t, :) = cmua(n).s(t).wave(:)';%, cmua(t).prm(4:6)];
        end
        [~, score] = pca(err);
        sPrmNrm = vertcat(cmua(n).s(:).prm);
        sPrm = [score(:, 1:6), sPrmNrm(:, 6)];
        
%         %2) simple parameters
%         sPrm = vertcat(cmua(n).s(:).prm);%spikes of current channel
        
        sPrmNrm = sPrm;%normalization
        for t = 1:size(sPrm, 2) %run over parameters
            sPrmNrm(:, t) = sPrm(:, t) - min(sPrm(:, t));
            sPrmNrm(:, t) = sPrmNrm(:, t) / max(sPrmNrm(:, t));
        end

        klustDir = 'c:\Users\Public\Distr\Утилиты_Редакторы\cluster\';%KlustaKwik direcotry
        klustKw = 'KlustaKwik-1_7.exe';%program filename
        fN = [klustDir, '_c', num2str(n)];%base name for temporary file

        maxClustN = '25';%initial number of possible clusters
        maxPossClustN = num2str(str2double(maxClustN) + 10);%maximal possible cluster number
        nStart = num2str(str2double(maxClustN) - 8);%start clusters number
        features = ones(size(sPrm, 2), 1);
        dlmwrite([fN '.fet.1'], size(sPrm, 2), 'delimiter', ' ', 'newline', 'pc')
        dlmwrite([fN '.fet.1'], sPrmNrm, 'delimiter', ' ', 'newline', 'pc', '-append')
        system([klustDir, klustKw, ' ', fN, ' 1', ...
            [' -MinClusters 0  -MaxClusters ', maxClustN,' -nStarts ', nStart, ...
            ' -MaxPossibleClusters ', maxPossClustN, ' -Screen 0'], ...
            ' -UseFeatures ' num2str(features)']);%, '-echo');
        idx = dlmread([fN '.clu.1']);
        idx(1) = [];%delete number of clusters
        idx = idx - 1;%cluster №1 - noise (become №0)

        disp(['===== ch', num2str(cmua(n).ch), ' ====='])
        disp([num2str(spkN), ' spikes'])
        disp([num2str(numel(unique(idx)) - 1), ' clusters'])

        for t = 1:spkN %run over spikes of current channel
            cmua(n).s(t).clst = idx(t);%cluster number
        end
    else
        disp(['===== ch', num2str(cmua(n).ch), ' ====='])
        disp('no spikes')
        disp('no clusters')
    end
end
%matlabpool close;

%% crosscluster histogram
clf
p1 = [10, 22];%[channel1, cluster1]
p2 = [10, chcl(10).c(11)];%[channel2, cluster2]
clstmua = [vertcat(cmua(p1(1)).s(:).tS) * zavp.siS * 1e3, vertcat(cmua(p1(1)).s(:).clst)];%[timestamp, cluster]
spk1 = clstmua((clstmua(:, 2) == p1(2)), 1);%spikes of cluster1
kk(1) = numel(unique(clstmua(:, 2))) - 1;

if all(p1 == p2) %autocorrelation
    spk2 = spk1;
    kk(2) = kk(1);
    bins = -10:1:10;%ms
    subplot(2, 1, 1)
else %crosscorrelation
    clstmua = [vertcat(cmua(p2(1)).s(:).tS) * zavp.siS * 1e3, vertcat(cmua(p2(1)).s(:).clst)];%[timestamp, cluster]
    spk2 = clstmua((clstmua(:, 2) == p2(2)), 1);%spikes of cluster2
    kk(2) = numel(unique(clstmua(:, 2))) - 1;
    bins = -20:1:20;%ms
    
%     %union of single-unit clusters
%     spk1 = [spk1; spk2];
%     spk2 = spk1;
%     bins = -10:1:10;%ms
end

spkDens = zeros(numel(bins), numel(spk1));
for t = 1:size(spk1, 1)
    spk0 = spk2(:, 1) - spk1(t, 1);
    spk0((spk0 > -0.5) & (spk0 < 0.5)) = [];%delete central event
    spkDens(:, t) = histc(spk0, bins);
end
bar(bins(1:(end - 1)) + diff(bins) / 2, sum(spkDens(1:(end - 1), :), 2), 'k')
xlim(bins([1, end]))
if all(p1 == p2)
    title(['ch=', num2str(cmua(p1(1)).ch), ' k=', num2str(p1(2)), '(', num2str(kk(1)), ')', '; n=', num2str(numel(spk1))])
else
    title(['ch', num2str(cmua(p1(1)).ch), ' k', num2str(p1(2)), '(', num2str(kk(1)), ')', '; n1=', num2str(numel(spk1)), ...
        ' | ch', num2str(cmua(p2(1)).ch), ' k', num2str(p2(2)), '(', num2str(kk(2)), ')', '; n2=', num2str(numel(spk2))])
end
% ylim([0, 0.005])

binStep = [3, 5];%[bin size (ms), precision (must be integer)]
%correlation with stimulus
if all(p1 == p2) %autocorrelation
    subplot(2, 1, 2), hold on
    bins = -100:5:200;%ms (stimulus-triggered)
    
    spk1 = (zavp.realStim(1).r - 1) / zavp.rarStep;%stimulus (ON-stim)
    spkDens = [];%zeros(numel(bins), numel(spk1));
    for t = 1:size(spk1, 1)
        spk0 = spk2(:, 1) - spk1(t, 1);
        spk0((spk0 > -0.5) & (spk0 < 0.5)) = [];%delete central event
        spkDens = [spkDens; spk0((spk0 >= bins(1)) & (spk0 <= bins(end)))];%(:, t) = histc(spk0, bins);
    end
    [s, basbins] = ZavOverHist(spkDens, binStep, bins([1, end]));
    plot(basbins, s, 'k.-')
    %plot(bins(1:(end - 1)) + diff(bins) / 2, sum(spkDens(1:(end - 1), :), 2), 'k.-')
    xlim(bins([1, end]))
	title(['ch=', num2str(cmua(p1(1)).ch), ' k=', num2str(p1(2)), '(', num2str(kk(1)), ')', '; n=', num2str(numel(spk1))])
    
    spk1 = (zavp.realStim(1).f - 1) / zavp.rarStep;%stimulus (OFF-stim)
    spkDens = [];%zeros(numel(bins), numel(spk1));
    for t = 1:size(spk1, 1)
        spk0 = spk2(:, 1) - spk1(t, 1);
        spk0((spk0 > -0.5) & (spk0 < 0.5)) = [];%delete central event
        spkDens = [spkDens; spk0((spk0 >= bins(1)) & (spk0 <= bins(end)))];%(:, t) = histc(spk0, bins);
    end
    [s, basbins] = ZavOverHist(spkDens, binStep, bins([1, end]));
    plot(basbins, s, '.-', 'Color', 0.5 * ones(1, 3))
    %plot(bins(1:(end - 1)) + diff(bins) / 2, sum(spkDens(1:(end - 1), :), 2), '.-', 'Color', 0.5 * ones(1, 3))
    xlim(bins([1, end]))
	title(['ch=', num2str(cmua(p1(1)).ch), ' k=', num2str(p1(2)), '(', num2str(kk(1)), ')', '; ', num2str(numel(spk1)), 'stim'])
end

%% plot multichannel clusters
clf
p1 = lL + rL + 1;
%p2 = [min(min(min(cat(3, cmua(:).wave)))), max(max(max(cat(3, cmua(:).wave))))];
chcl = [[10, 12]; [10, 22]; [11, 3]];%channel-cluster
for m = 1:size(chcl, 1) %run over clusters
    n = chcl(m, 1);%channel(index)
    k = chcl(m, 2);%cluster
    clstmua = [vertcat(cmua(n).s(:).tS) * zavp.siS * 1e3, vertcat(cmua(n).s(:).clst)];%[timestamp, cluster]
    jj = find((clstmua(:, 2) == k));
    spkN = numel(jj);%number of requested spikes
    spkDens = zeros(p1, spkN);
    ii = 1:(floor(log2(spkN) * 1.5) + 1):spkN;%random spikes
    for z = 1:3 %run over channels
        for t = 1:spkN %run over spikes
            spkDens(:, t) = cmua(n).s(jj(t)).wave(:, z);
        end
        subplot(3, size(chcl, 1), m + ((z - 1) * size(chcl, 1))), hold on
        if (z == 1)
            title(['ch', num2str(cmua(n).ch), ' k=', num2str(k), '; n=', num2str(spkN)])
        end
        plot(-lL:rL, spkDens(:, ii), 'Color', 0.8 * ones(1, 3)), xlim([-lL, rL])
        plot(-lL:rL, mean(spkDens, 2), 'Color', zeros(1, 3), 'LineWidth', 3)
        ylim([-90, 30])
        %ylim([min(mean(spkDens, 2)), max(mean(spkDens, 2))])
    end
    %disp([diff(ii(1:2)), length(ii), spkN])
end



%% === spikes distribution analysis === %%
clf
ch = 8;%scanning channel
z = 2;%type of spikes distribuion analysis

switch z %spikes frequency
case 1
    spkNew(1:hd.lActualEpisodes) = struct('s', []);
    spkDens = [];%pulse-to-pulse intervals
    for sw = 1:hd.lActualEpisodes
        spkNew(sw).s = spks(ch, sw).tStamp * 1e-3;%spikes moments (s)
        spkDens = [spkDens; diff(spks(ch).tStamp(spks(ch).tStamp <= tW))];%pulse-to-pulse intervals
    end
    hf = 0:5000;%frequency range
    hs = hist(1e3 ./ spkDens, hf);
    inds = find(hs(1:end - 1) > 0);
    subplot(2, 1, 1)
    plot(hf(inds), hs(inds), '.-')
    subplot(2, 1, 2)
    [spkCnt, bins, err] = isi(spkNew, [1, 1.2], 1, 100);%chronux            
    %plot(bins, spkCnt)
    errorbar(bins, spkCnt, err, err)
case 2 %sparcification index
    k = 1;%number of interval of spontaneous activity
    overlapCoef = zeros(hd.lActualEpisodes * (numel(zavp.realStim(1).r) + 1), 1);%overlapping coefficients
    for sw = 1:hd.lActualEpisodes
        if (numel(zavp.realStim(sw).r) <= 0)%no stimulus
            p2 = size(lfp, 1);%length of spontaneous activity
            spkDens = diff(spks(ch, sw).tStamp(spks(ch, sw).tStamp <= p2));
            tW = p2 / (numel(spkDens) + 1);%average time for single spike
            for jj = 1:numel(spkDens)
                overlapCoef(k) = overlapCoef(k) + min(tW, spkDens(jj));
            end
            overlapCoef(k) = overlapCoef(k) / p2;
            k = k + 1;%number of interval of spontaneous activity
        else%stimulus
            p1 = 0;%begin of record (ms)
            for t = 1:numel(zavp.realStim(sw).r)
                p2 = (zavp.realStim(sw).r(t) / zavp.rarStep);%moment of next stimulus (ms)
                if ((p2 - p1) > 900)%interval of silence quite long
                    spkDens = diff(spks(ch, sw).tStamp((spks(ch, sw).tStamp >= p1) & (spks(ch, sw).tStamp <= p2)));
                    tW = (p2 - p1) / (numel(spkDens) + 1);%average time for single spike
                    for jj = 1:numel(spkDens)
                        overlapCoef(k) = overlapCoef(k) + min(tW, spkDens(jj));
                    end
                    overlapCoef(k) = overlapCoef(k) / (p2 - p1);
                    k = k + 1;%number of interval of spontaneous activity
                end
                p1 = p2 + 1e3;%teoretical begin of next silence interval
            end
        end
    end
    overlapCoef(k:end) = [];%delete empty
    disp(['overlapCoef = ', num2str(mean(overlapCoef))]);
    
    k = 1;%number of interval of spontaneous activity
    spkPerSec = zeros(hd.lActualEpisodes * (numel(zavp.realStim(1).r) + 1), 1);%spkes per second
    for sw = 1:hd.lActualEpisodes
        if (numel(zavp.realStim(sw).r) <= 0)%no stimulus
            p2 = size(lfp, 1);%length of spontaneous activity
            spkPerSec(k) = sum(spks(ch, sw).tStamp <= p2) / (p2 * 1e-3);
            k = k + 1;%number of interval of spontaneous activity
        else
            p1 = 0;%begin of record (ms)
            for t = 1:numel(zavp.realStim(sw).r)
                p2 = (zavp.realStim(sw).r(t) / zavp.rarStep);%moment of next stimulus (ms)
                if ((p2 - p1) > 900)%interval of silence quite long
                    spkPerSec(k) = sum((spks(ch, sw).tStamp >= p1) & (spks(ch, sw).tStamp <= p2)) / ((p2 - p1) * 1e-3);
                    k = k + 1;%number of interval of spontaneous activity
                end
                p1 = p2 + 1e3;%teoretical begin of next silence interval
            end
        end
    end
    spkPerSec(k:end) = [];%delete empty
    disp(['spkPerSec = ', num2str(mean(spkPerSec))]);
case 3 %power spectrum for spikes
    prm = zavp.prm; prm.Fs = 1 / zavp.siS;%original discretization frequency
    prm.fpass = [0, 200];
    prm.trialave = 1;
    [spkSpct, spkF] = mtspectrumpt(allSpk(1:100) * 1e-3, prm);
    plot(spkF, spkSpct);
end

        
%% === man EEG correlation in time === %%
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [-10e3, 10e3];%left and right shifts from synchro-point (ms)
rCh = [3, 5];%1:hd.nADCNumChannels;%channels to be read and draw
sn = 12;%1:size(segms, 1);%wanted sweeps or segments

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)

if rawData, prm.Fs = 1 / zavp.siS; ...%sampling frequency for raw data
else prm.Fs = 1 / (zavp.siS * zavp.rarStep); end%sampling frequency for rare data
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)

tW = 1e3;%moving window (ms)
corCf = zeros(size(whtLFP, 1), 1);%correlation coefficient
for t = 1:(size(whtLFP, 1) - tW)
    chr = corrcoef(whtLFP(t:(t + tW), :));
    corCf(t) = chr(1, 2);%
end
clf
plot(tm, [corCf, whtLFP / max(max(whtLFP)) - 1])


%% === circular statistics === %%
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [-500, 500];%left and right shifts from synchro-point (ms)
rCh = 2;%channels for LFP
sn = 1:size(segms, 1);%wanted sweeps or segments

whtLFP = ZavFilter(lfp, prm.Fs, 'high', 2, 2);%filtered LFP
whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, whtLFP, rCh, rawData, nlxVer);%lfp phased with respect to synchro-events
whtLFP = squeeze(whtLFP);%average by sweeps (segments)

k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)

ch = rCh;%channel for spikes

sw = 1;
clear spkNew
spkNew(1:size(whtLFP, 2)) = struct('s', [], 'ph', []);%spikes phases
k = find(tm >= 0, 1, 'first');
lfpPh = zeros(size(whtLFP));
for t = 1:size(whtLFP, 2)
    lfpPh(:, t) = unwrap(angle(hilbert(-smooth(whtLFP(:, t), 21))));%hilbert phase
%     p1 = find(diff(lfpPh(k:-1:1, t)) > pi, 1, 'first'); p1 = k - p1;
%     p2 = find(diff(lfpPh(k:end, t)) < -pi, 1, 'first'); p2 = k + p2;
%     lfpPh([1:p1, p2:end], t) = Inf;
    lfpPh(:, t) = lfpPh(:, t) - lfpPh(k, t);
    lfpPh(lfpPh(:, t) < (-4 * pi), t) = Inf;
    lfpPh(lfpPh(:, t) > (4 * pi), t) = Inf;
    tW = segms(t, 1) + segmEdge;
    spkNew(t).s = spks(ch, sw).tStamp(spks(ch, sw).tStamp > tW(1) & spks(ch, sw).tStamp < tW(2)) - tW(1);
    spkNew(t).ph = zeros(size(spkNew(t).s));
    z = size(spkNew(t).s, 1);
    for g = 1:z
        p1 = max(1, floor(spkNew(t).s(g)));
        p2 = min(ceil(spkNew(t).s(g)), size(tm, 2));
        if p1 ~= p2
            spkNew(t).ph(g) = interp1([p1, p2], lfpPh([p1, p2], t), spkNew(t).s(g), 'linear');%phase of spike
        else
            spkNew(t).ph(g) = lfpPh(p1, t);%phase of spike
        end
    end
end
clf
p1 = vertcat(spkNew(:).ph); p1 = p1(isfinite(p1)); hist(p1, (-4 * pi):0.1:(4 * pi))


%% === parameters in time (evoked) === %%
ch = 8;
spkDens = zeros(size(segms, 1), 2);
sw = 1;
for t = 1:size(segms, 1)
    if (isfinite(sepOnsetPeak(t, 1)) && (sepOnsetPeak(t, 1) > 0))
        segmEdge = segms(t, 1) + sepOnsetPeak(t, 1) + [0, 20];%SEP period
        spkDens(t, 1) = sum((spks(ch, sw).tStamp >= segmEdge(1)) & (spks(ch, sw).tStamp <= segmEdge(2)));%mua during SEP
        segmEdge = segms(t, 1) + sepOnsetPeak(t, 1) + [20, 520];%period after SEP (burst)
        spkDens(t, 2) = sum((spks(ch, sw).tStamp >= segmEdge(1)) & (spks(ch, sw).tStamp <= segmEdge(2)));%mua after SEP (burst)
    end
end

binStep = [240, 120];%[bin size, time step] (sec)
tm = hd.inTTL_timestamps.t(:, 1) * 1e-6;%actual time (sec)
bins = 0:binStep(2):diff(hd.recTime);%requested moments
p = zeros(numel(bins), 5);%parameter values
for t = 1:numel(bins)
    jj = (tm > bins(t)) & (tm <= (bins(t) + binStep(1)));
    p(t, 1:3) = mean(sepOnsetPeak(jj, :), 1);%SEP
    p(t, 4:5) = mean(spkDens(jj, :), 1);%MUA
end
clf
subplot(3, 1, 1)
plot(bins, p(:, 3), '.-')%amplitude (time in seconds)
subplot(3, 1, 2)
plot(bins, p(:, 1:2), '.-')%latency
subplot(3, 1, 3)
plot(bins, p(:, 4:5), '.-')%MUA during and after SEP


%% === parameters in time (spontaneous) === %%
binStep = [240, 120];%[bin size, time step] (sec)
bins = 0:binStep(2):diff(hd.recTime);%requested moments
p = Inf(numel(bins), 5);%parameter values
if exist([zavp.file(1:(end - 1)), '_brst.mat'], 'file')
    load([zavp.file(1:(end - 1)), '_brst.mat'])
    tm = brst(:, 1) * 1e-3;%actual time (sec)
    for t = 1:numel(bins)
        jj = (tm > bins(t)) & (tm <= (bins(t) + binStep(1)));
        p1 = brst(jj, :);%bursts in segment
        if (size(p1, 1) > 0)
            p1 = p1(p1(:, 3) >= 5, :);%bursts with MUA
            p(t, 1) = size(p1, 1) / (binStep(1) / 60);%bursts with MUA (1 / min)
            p(t, 2) = mean(round(1e4 * p1(:, 3) ./ (diff(p1(:, 1:2), 1, 2) * 1e-3)) * 1e-4);%MUA in bursts (1 / s)
            p(t, 3) = sum(jj) / (binStep(1) / 60);%total bursts (1 / min)
        end
    end
end

if exist([zavp.file(1:(end - 1)), '_mbs.mat'], 'file')
    load([zavp.file(1:(end - 1)), '_mbs.mat'])
    tm = mbs(:, 1) * 1e-3;%actual time (sec)
    for t = 1:numel(bins)
        jj = (tm > bins(t)) & (tm <= (bins(t) + binStep(1)));
        if (sum(jj) > 0)
            p(t, 4) = sum(jj) / (binStep(1) / 60);%movements (1 / min)
        end
    end
end
ch = 8;
tm = spks(ch).tStamp * 1e-3;%actual time (sec)
for t = 1:numel(bins)
    jj = (tm > bins(t)) & (tm <= (bins(t) + binStep(1)));
    p(t, 5) = sum(jj) / binStep(1);%MUA / s
end

clf
subplot(3, 1, 1)
plot(bins, p(:, [1, 3]), '.-')%bursts frequency
subplot(3, 1, 2)
plot(bins, p(:, 2), '.-')%MUA in bursts (1 / 2)
subplot(3, 1, 3)
plot(bins, p(:, 5), '.-')%movements frequency

%% === sort troughs to classes === %%
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [-200, 200];%left and right shifts from synchro-point (ms)
rCh = 8 + (-5:4);%channels to be read and draw
sn = 1:10;%1:size(segms, 1);%wanted segments (practically for detected troughs)
%sn = 1;%1:hd.lActualEpisodes;%sweeps to be read

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
%whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfpEx, rCh, rawData, nlxVer)
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)

prm.fpass = [20 100];
prm.tapers = [2, 3];%tapers
prm.pad = 1;%padding factor
prm.trialave = 1;

trCls = 20;%number of troughs classes

aGrp = zeros(trCls, 2);%groups of troughs amplitudes
clear tGP
tGP(1:trCls) = struct('z1', 0, 'z2', 0, 'lfp', zeros(201, 100), ...
    'chr', struct([]), 'f', [], 'phi', struct([]), 'tSh', [], 'r', struct([]));%parameters for groups of troughs
jj = [];
for ch = 8 %run over channels
    for sw = 1:hd.lActualEpisodes
        jj = [jj; respPrm(ch, sw).t(:, 4)];%all amplitudes
    end
end
%p2 = floor(length(jj) / trCls);%number of trough in each group
p2 = (max(jj) - min(jj)) / trCls;%amplitude step

aGrp(1, 1:2) = min(jj);%left-right board of amplitude range
for t = 1:trCls %run over troughs groups
    if (t > 1)
        aGrp(t, 1:2) = aGrp(t - 1, 2);%left-right board of amplitude range
    end
    %aGrp(t, 2) = aGrp(t, 1) + p2;
    while (sum((jj >= aGrp(t, 1)) & (jj < aGrp(t, 2))) < p2)
        aGrp(t, 2) = aGrp(t, 2) + 0.1;%right board of amplitude range
        if (aGrp(t, 2) >= 0)
            break;%out from while (sum((brst(ch...
        end
    end
    tGP(t).tSh = 100:-1:-100;%time shifts range
    for c = 3:(hd.nADCNumChannels - 2)
        tGP(t).r(c).r = zeros(length(tGP(t).tSh), 100);
        tGP(t).chr(c).chr = [];
        tGP(t).phi(c).phi = [];
    end
    tGP(t).spk(1:hd.nADCNumChannels, 1:100) = struct('s', []);
    tGP(t).s(1:hd.nADCNumChannels) = struct('s', []);
end

for ch = 8 %run over channels
    for sw = 1:hd.lActualEpisodes
        for p2 = 1:size(respPrm(ch, sw).t, 1) %run over troughs
            t = find(((respPrm(ch, sw).t(p2, 4) >= aGrp(:, 1)) & (respPrm(ch, sw).t(p2, 4) < aGrp(:, 2))) == 1, 1, 'first');
            if ~isempty(t)
                tW = respPrm(ch, sw).t(p2, 2) + segms([], 1) + [-100, 100];
                if ((tW(1) > 0) && (tW(2) < size(whtLFP, 1)))
                    tGP(t).z1 = tGP(t).z1 + 1;%%numbers of troughs in the group
                    tGP(t).lfp(1:201, tGP(t).z1) = whtLFP(tW(1):tW(2), ch, sw);%lfp segment
                    for c = 3:(hd.nADCNumChannels - 2)
                        tGP(t).spk(c, tGP(t).z1).s = (spks(c, sw).s(spks(c, sw).s >= tW(1) & spks(c, sw).s <= tW(2)) - tW(1)) / 1e3;
                        tGP(t).s(c).s = [tGP(t).s(c).s; tGP(t).spk(c, tGP(t).z1).s];%union of all spikes
                    end

                    if (~isempty(tGP(t).lfp) && ~isempty(tGP(t).spk) && ((tW(1) + tGP(t).tSh(1)) > 0) && ((tW(2) + tGP(t).tSh(end)) < size(whtLFP, 1)))
                        tGP(t).z2 = tGP(t).z2 + 1;
                        for c = 5:(hd.nADCNumChannels - 2)
                            [~, cS] = ZavPointToContin(0:5e-3:0.2, hist(tGP(t).spk(c, tGP(t).z1).s, 0:5e-3:0.2), 1e-3);
                            for z = 1:length(tGP(t).tSh)
                                if (((tW(1) + tGP(t).tSh(z)) > 0) && ((tW(2) + tGP(t).tSh(z)) < size(whtLFP, 1)))
                                    p1 = corrcoef(cS, whtLFP((tW(1):tW(2)) + tGP(t).tSh(z), ch, sw));
                                else
                                    p1 = [0, 0];
                                end
                                if any(isnan(p1))
                                    p1 = [0, 0];
                                end
                                tGP(t).r(c).r(z, tGP(t).z2) = p1(1, 2);%correlation coefficient
                            end
                        end
                    end
                end
            end
        end
    end
end

for t = 1:trCls %run over troughs groups
    tGP(t).lfp(:, tGP(t).z1 + 1:end) = [];%delete excess
    tGP(t).spk(:, tGP(t).z1 + 1:end) = [];%delete excess
    for ch = 3:(hd.nADCNumChannels - 2)
        tGP(t).r(ch).r(:, (tGP(t).z2 + 1):end) = [];%delete excess
        if (~isempty(tGP(t).lfp) && ~isempty(tGP(t).spk(ch, :)))
            [chr, phi, ~, ~, ~, f] = coherencycpt(tGP(t).lfp, tGP(t).spk(ch, :), zavp.prm);
            tGP(t).f = f; tGP(t).chr(ch).chr = chr; tGP(t).phi(ch).phi = phi;
        end
    end
end

ch = 5;%scanning channel
bins = 0:2e-3:0.2;
psp = zeros(length(bins), trCls);%spikes density class by class
mLCrl = zeros(length(tGP(1).tSh), trCls);%spikes-lfp cross-correlation class by class
for t = 1:trCls %run over troughs groups
    psp(:, t) = hist(tGP(t).s(ch).s, bins) / (size(tGP(t).lfp, 2) + 0);
    mLCrl(:, t) = mean(tGP(t).r(ch).r, 2);
end
clf
pcolor(bins - 0.1, 1:trCls, psp')
%pcolor(-100:100, 1:trCls, abs(mLCrl'))
shading flat; colormap('default'); colorbar; hold on
plot([0, 0], [1, trCls], 'r', 'LineWidth', 2)
% fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
% title(fN)

%% === lfp-movement correlation === %%
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [-200, 200];%left and right shifts from synchro-point (ms)
rCh = 8 + (-5:4);%channels to be read and draw
sn = 1:10;%1:size(segms, 1);%wanted segments (practically for detected troughs)
%sn = 1;%1:hd.lActualEpisodes;%sweeps to be read

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
%whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfpEx, rCh, rawData, nlxVer)
whtLFP = mean(whtLFP, 3);%average by sweeps (segments)
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)

lenThrsh = 250;%threshold of movement length
znak = 1;%sign (1 - larger than lenThrsh; -1 - smaller than lenThrsh)
z = 2;%1 - correlation (spikes-movement or lfp-movement)
      %2 - start-movement - start-lfp interval histogram
if (z == 1)%correlation (spikes-movement or lfp-movement) %%%
    binStep = 100;%bin step (ms)
    bins = (binStep / 2):binStep:(p1(2) - p1(1) + 1);%centrs of interval (binning)
    jj = -300:300;%time shifts range
    r = zeros(length(jj), size(mbs, 1));
    k = 1;%counter of good segments
    for t = 1:size(mbs, 1)
        tW = mbs(t, 1) + p1;%start-stop time (seconds)
        %spkNew = spks(ch).s((spks(ch).s >= tW(1)) & (spks(ch).s <= tW(2))) - tW(1);%spikes time (seconds)
        %spkS = smooth(interp1(bins, hist(spkNew, bins), 1:(p1(2) - p1(1) + 1)), 27);
        if (((znak * mbs(t, 4)) > (znak * lenThrsh)) && ...
            ((tW(1) + jj(1)) > 0) && ((tW(2) + jj(end)) < size(whtLFP, 1)))
            for z = 1:length(jj)
                p2 = corrcoef(whtLFP(tW(1):tW(2), 1), whtLFP((tW(1):tW(2)) + jj(z), 2));
                %p2 = corrcoef(spkS, whtLFP((tW(1):tW(2)) + jj(z), 2));
                if any(isnan(p2))
                    p2 = [0, 0];
                end
                r(z, t) = p2(1, 2);%correlation coefficient
            end
            k = k + 1;%counter of good segments
        end
    end
    r(:, k:end) = [];%delete empty
    clf
    plot(jj, mean(r, 2), 'k', 'LineWidth', 1)
elseif (z == 2)%start-movement - start-lfp interval histogram %%%
    bins = -200:20:500;
    k = 1;%counter of good segments
    p2 = zeros(size(segms, 1), 1);
    for t = 1:size(mbs, 1)
        tW = mbs(t, 1) + p1;%start-stop time (seconds)
        if (((znak * mbs(t, 4)) > (znak * lenThrsh)) && ...
            (tW(1) > 0) && (tW(2) < size(whtLFP, 1)))
            %jj = find(whtLFP(tW(1):tW(2)) < r, 1, 'first');
            jj = brst(ch).t(2).t((brst(ch).t(2).t(:, 2) > tW(1)) & (brst(ch).t(2).t(:, 2) < tW(2)), 2);
            if ~isempty(jj)
                %p2(k) = jj + tW(1) - 1 - mbs(t, 1);
                p2(k) = jj(1) - mbs(t, 1);
                k = k + 1;%counter of good segments
            end
        end
    end
    p2(k:end) = [];%delete empty
    clf
    hist(p2, bins)
end

% fN = regexprep(fN, '_', '\\_');%replace character '_' by '\_'
% title(fN)

%% === fluctuations === %%
%find changes of parameters of naturaly randome processes
rawData = false;%if raw data needed then rawData = 1(true)

segmEdge = [0, 300];%left and right shifts from synchro-point (ms)
rCh = 8;%1:hd.nADCNumChannels;%channels to be read and draw

t = 1;
ajj = (0:100:(size(lfp, 1) - diff(segmEdge)))';
segms = zeros(numel(ajj) * hd.lActualEpisodes, 4);
for sw = 1:hd.lActualEpisodes
    jj = t:(t + numel(ajj) - 1);
    segms(jj, 1) = ajj;%position of synchro-event (ms)
    segms(jj, 2) = rCh;%number of channel where syncro-event was detected
    segms(jj, 3) = sw;%number of sweeps where syncro-event was detected
    segms(jj, 4) = -1;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
    t = t + numel(ajj);
end
sn = 1:size(segms, 1);%wanted sweeps or segments

whtLFP = ZavSynchLFP(zavp, hd, segms(sn, :), segmEdge, lfp, rCh, rawData, nlxVer);%lfp phased with respect to stimulus moments
%whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfpEx, rCh, rawData, nlxVer)

whtLFP = squeeze(whtLFP);%average by sweeps (segments)
k = 1 + ((zavp.rarStep - 1) * double(rawData == true));%multiplier for time matrix
tm = ((segmEdge(1) * k):(segmEdge(2) * k)) / k;%time matrix for raw or downsampled data (ms)

% clf, hold on
bins = min(lfp(:, rCh)):10:(max(lfp(:, rCh)));
ajj = zeros(numel(sn), length(bins));
for t = 1:numel(sn)
    ajj(t, :) = smooth(hist(whtLFP(:, t), bins), 8);
end

%% === channel to channels spikes histogram === %%
rCh = [7; 11];%first row - first group, second row - second group
spk1 = vertcat(spks(rCh(1, :)).s);%first spikes group
spk1 = sort(spk1);
spk2 = vertcat(spks(rCh(2, :)).s);%first spikes group
spk2 = sort(spk2);

binStep = 1;%bin step (ms)
bins = 0:binStep:500;%binning of time segment

jj = (spk2 > (spk1(1) + bins(1))) & (spk2 < (spk1(1) + bins(end)));
spkHist = hist(spk2(jj) - spk1(1), bins);
for t = 2:numel(spk1)
    jj = (spk2 > (spk1(t) + bins(1))) & (spk2 < (spk1(t) + bins(end)));
    spkHist = spkHist + hist(spk2(jj) - spk1(t), bins);
end
spkHist = spkHist / numel(spk1);
%spkHist([1, end]) = 0;

clf
plot(bins, spkHist, 'k')



otherwise
    disp('unexpected operation mode')
end

%% === StimMoments (function) === %%
function [stimCh, realStim] = StimMoments(flNm, hd, isStim, synchro)
%[stimCh, realStim] = StimMoments(flNm, hd, isStim)
%define moments of stimulus
%
%INPUTS
%flNm - full pathname of file
%hd - header of data file
%isStim - data contain stimulus artefact
%synchro - known moment of synchro-impulses
%
%OUTPUS
%stimCh - number of channel containing synchro-impulses
%realStim - moments of detected synchro-impulses (samples (from 1))

stimCh = NaN;%no stimulus channel initially
%rRealStim(1:hd.lActualEpisodes) = struct('r', [], 'f', []);%rare version of zavp.realStim (ms)
%rRealStim(sw).r|f = round(zavp.realStim(sw).r|f / zavp.rarStep);%rare version of zavp.realStim.r|f (ms)

if isStim %try to find synchro impulses
    if (isequal(hd.fFileSignature, 'Neuralynx') && ... neuralynx file
            ~isempty(vertcat(synchro(:).t)))%synchro-array exist
        realStim(1:length(synchro)) = struct('r', [], 'f', []);%stimulation moments (samples, rise and fall fronts (from 1))
        for ch = 1:length(synchro)
            realStim(ch).r = synchro(ch).t(:, 1);%rise front of synchro-impulse (samples from sweep start (from 1))
            realStim(ch).f = synchro(ch).t(:, 2);%fall front of synchro-impulse (samples from sweep start (from 1))
        end
    elseif (isequal(hd.fFileSignature, 'eeg EDF+') && ~isempty(vertcat(synchro(:).t)))%man EEG
        realStim(1).r = zeros(size(synchro, 1), 1);%rise front of synchro-impulse (samples from sweep start (from 1))
        realStim(1).f = zeros(size(synchro, 1), 1);%fall front of synchro-impulse (samples from sweep start (from 1))
        for sw = 1:size(synchro, 1)
            realStim(1).r(sw) = synchro(sw, 1);%stimulation moments (rise front, samples from sweep start (from 1))
            realStim(1).f(sw) = synchro(sw, 2);%ends of stimulus (fall front, samples from sweep start (from 1))
        end
    elseif (isequal(hd.fFileSignature(1:3), 'ABF') || isequal(hd.fFileSignature, 'DAQ')) %ABF or DAQ
%         if isequal(hd.fFileSignature(1:3), 'ABF') %ABF
%             data = ZavLoadData(flNm, hd, 'a', 1, 0, 'e', nlxVer);%min(120, diff(hd.recTime)));%one sweep (segment) for stimulus channel definition
%         else %DAQ
%             if ((hd.ObjInfo.SamplesPerTrigger == Inf) || isequal(hd.ObjInfo.SamplesPerTrigger, 'Inf')) %gap free
%                 data = ZavLoadData(flNm, hd, 'a', min(2, hd.lActualEpisodes), 0, min(20, diff(hd.recTime)), nlxVer);%one segment for stimulus channel definition
%             else %sweeps
%                 data = ZavLoadData(flNm, hd, 'a', min(2, hd.lActualEpisodes), 0, min(20, 2 * diff(hd.recTime) / hd.lActualEpisodes), nlxVer);%one sweep for stimulus channel definition
%             end
% %             difData = diff(data, 1, 1);%differential
%         end
%         rws = Inf(4, hd.nADCNumChannels); rws(2, :) = -Inf;
        mdfp = 1200 / hd.si;%window (points corresponding to 1.2 ms) minimal duration of stimulus
%         for ch = 1:hd.nADCNumChannels %[1:3, (hd.nADCNumChannels - 2):hd.nADCNumChannels] %run over channels
%             for t = 1:round(mdfp / 5):(size(data, 1) - mdfp)
%                 p = max(abs(data(t:(t + mdfp), ch))) - min(abs(data(t:(t + mdfp), ch)));%difference between maximum and minimum
% %                 p(2) = max(abs(difData(t:(t + mdfp), ch))) - min(abs(difData(t:(t + mdfp), ch)));%difference between maximum and minimum
%                 if (p < 1e-9)
%                     p = 1e-6;
%                 end
%                 if (p < rws(1, ch))
%                     rws(1, ch) = p;%minimal difference
%                 end
%                 if (p > rws(2, ch))
%                     rws(2, ch) = p;%maximal difference
%                     rws(3, ch) = t;%point
%                 end
% %                 if (p(2) < rws(3, ch))
% %                     rws(3, ch) = p(2);%minimal difference
% %                 end
% %                 if (p(2) > rws(4, ch))
% %                     rws(4, ch) = p(2);%maximal difference
% %                 end
%             end
%             rws(5, ch) = ch;%number of channel
%         end
%         rws(4, :) = rws(2, :) ./ rws(1, :);%ration of minimal and maximal difference
%         rws = rws(:, isfinite(rws(4, :))); [~, t] = max(rws(4, :));
% 
%         stimCh = rws(5, t);%channel with stimulus artefact 
        stimCh = 17;%special for daq of Marseille_2017

        realStim(1:hd.lActualEpisodes) = struct('r', [], 'f', []);%stimulation moments (samples, rise and fall fronts (from 1))
        data = squeeze(ZavLoadData(flNm, hd, stimCh, 'a', 0, 'e', []));%load stimulus markers
        for sw = 1:hd.lActualEpisodes %run over sweeps
            stmThrsh = max(abs(data(:, sw))) / 3;%threshold for detection of stimulus
            onStim = find(abs(data(:, sw)) > stmThrsh);%stimulation moments (samples from beginning of the sweep, rise front (from 1))
            offStim = onStim(end:-1:1);%end of stimulus (samples from beginning of the sweep, fall front (from 1))
            reps = find(diff(onStim) < mdfp);%repeat points (during stimulus)
            onStim(reps + 1) = [];%stimulation moments (samples from beginning of the sweep, rise front (from 1))
            reps = find(diff(offStim) > -mdfp);%repeat points (during stimulus)
            offStim(reps + 1) = [];%end of stimulus (samples from beginning of the sweep, fall front (from 1))
            offStim = offStim(end:-1:1);%end of stimulus (samples from beginning of the sweep, fall front (from 1))

            %common algorithm (but Marseille_2017)
            %realStim(sw).r = onStim;%stimulation moments (samples from beginning of the sweep, rise front (from 1))
            %realStim(sw).f = offStim;%end of stimulus (samples from beginning of the sweep, fall front (from 1))
            
            %special for daq of Marseille_2017
            realStim(sw).r = onStim(1:2:end);%stimulation moments (samples from beginning of the sweep, rise front (from 1))
            realStim(sw).f = onStim(2:2:end);%end of stimulus (samples from beginning of the sweep, fall front (from 1))
%             if ((realStim(sw).f(1) - realStim(sw).r(1)) > (2.2e5 / hd.si)) %wrong sequence
%                 offStim = realStim(sw).r;
%                 realStim(sw).r = realStim(sw).f;
%                 realStim(sw).f = offStim;
%             end
        end
    else
        realStim = struct('r', [], 'f', []);%empty
    end
else %no stimulation
    realStim = struct('r', [], 'f', []);%empty
end


