function [data, ttlIn, hd, spkTS, spkSM] = ZavNrlynx(pf, hd, rCh, nlxVer, strt, stp)
%[data, ttlIn, hd, spkTS, spkSM] = ZavNrlynx(pf, hd, rCh)%, strt, stp)
%read neuralynx. NeuralynxMatlabImportExport_v501 require
%
%INPUTS
%pf - pathname of directory with files *.ncs, *.nev, *.nse (or full pathename of one of file in the directory)
%hd - header (full, for all channels)
%rCh - numbers of channels to be read
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%strt - start time to read record (ms from beginning of record) - not used
%stp - stop time to read record (ms from beginning of record) - not used
%
%OUTPUTS
%data - signal samples (microvoltes)
%ttlIn - moments of synchro-TTL inputs (samples from sweep beginning)
%hd - file header (information about record)
%spkTS - spikes appearence moments (mks)
%spkSM - spikes samples

if (strcmp(pf((end - 3):end), '.ncs') || strcmp(pf((end - 3):end), '.nev') || strcmp(pf((end - 3):end), '.nse'))
    t = find(pf == '\', 1, 'last');%last slash
    pf = pf(1:t);%directory paht only
end
if (pf(end) ~= '\') %no slash
    pf(end + 1) = '\';
end
if ~exist(pf, 'dir') %&& ~exist(pf, 'file')
    [~, pf] = uigetfile('*.ncs; *.nev; *.nse', 'Select file', pf);%open dialog
end

rChNum = 0;%counter of channels
dirCnt = dir(pf);%directory content
ncsFiles(1:length(dirCnt)) = struct('f', '', 'bytes', -Inf, 'r', false);%full names of channel-files
%ii = true(1, length(dirCnt));
largestF = Inf;%number of file (in ncsFiles structure) with largest size
for t = 1:length(dirCnt) %run over files in NLX directory (ncs, we must scan all channels here)
    if ((~dirCnt(t).isdir) && (length(dirCnt(t).name) > 3)) %not directory and name good
        if isequal(dirCnt(t).name(end - 3:end), '.ncs') %data-file of cheetah
            recflnmlen = length(dirCnt(t).name);%length of filename
            z = recflnmlen - 4;%character position
            while (z > 1) %run over characters of filename
                c = double(dirCnt(t).name(z));%double representation of a character
                if ((c < 48) || (c > 57)) %not a number
                    break;%out of (while z)
                end
                z = z - 1;%next character
            end
            ch = str2double(dirCnt(t).name((z + 1):(recflnmlen - 4)));%number of channel (3-letteral name + number of channel)
            ncsFiles(ch).f = [pf, dirCnt(t).name];%full name of channel-file
            ncsFiles(ch).bytes = dirCnt(t).bytes;%size of file (bytes)
            %ii(ch) = false;
            if (~isfinite(largestF) || (ncsFiles(ch).bytes > ncsFiles(largestF).bytes))
                largestF = ch;%number of file (in ncsFiles structure) with largest size
            end
            rChNum = rChNum + 1;%increase counter of channels
        end
    end
end
ncsFiles(~isfinite(vertcat(ncsFiles(:).bytes))) = [];%delete empty
for t = 1:length(ncsFiles) %run over ncs-files
    if (ncsFiles(t).bytes == ncsFiles(largestF).bytes) %good file size
        ncsFiles(t).r = true;%we can read file fully
    end
end

if (isempty(rCh) || isequal(rCh, 'a')) %all channels requested
    if ~isempty(hd) %no header
        rChNum = hd.nADCNumChannels;
    end
    rCh = 1:rChNum;%read all channels
end
rChNum = length(rCh);%number of wanted channels

if isempty(hd) %no header
    cscHd = Nlx2MatCSC(ncsFiles(largestF).f, [0 0 0 0 0], 1, 1, []);%header (lfp)
    dsi = (1e6 / NlxParametr(cscHd, 'SamplingFrequency'));%sample interval (mks)
    cscTStmp = Nlx2MatCSC(ncsFiles(largestF).f, [1 0 0 0 0], 0, 1, []);%timestamps of recorded samples (read largest file)
else %we have a header
    dsi = hd.si;%sample interval (mks)
    if ~isfield(hd, 'cscTStmp') %no csc-timestamps saved
        cscTStmp = Nlx2MatCSC(ncsFiles(largestF).f, [1 0 0 0 0], 0, 1, []);%timestamps of recorded samples (read largest file)
    else %we have saved csc-timestamps
        cscTStmp = hd.cscTStmp;%timestamps of recorded samples (read largest file)
    end
end
if any(diff(cscTStmp, 2) > 5) %datasections lost
    disp('large error of timestamp')%datasections lost
end

%%% stimulus moments (if exist) %%%
fileToRead = [pf, 'Events.nev'];%pathname of file to be read
if exist(fileToRead, 'file') %events-file exist
    [evntTStmp, ttl, evntStr] = Nlx2MatEV(fileToRead, [1 0 1 0 1], 0, 1, []);%read events, timestamps and event strings
    mStStRec = FindStStStemp(cscTStmp, evntStr, evntTStmp, dsi);%find moments of "Starting Recording" and "Stopping Recording" (microseconds)
    if (nargout > 1) %synchro events requested
        inEvntOn = zeros(1, length(evntStr));%numbers of events when input ports changed
        for t = 1:length(evntStr) %run over event strings
            inEvntOn(t) = t * double(~isempty(strfind(evntStr{t}, 'Input')));%find input events
            %inEvntOn(t) = t * double(~isempty(strfind(evntStr{t}, 'Output')));%find output events
        end
        inEvntOff = inEvntOn((ttl <= 0) & (inEvntOn > 0));%number of events when input TTL ports are in state 'OFF'
        inEvntOn = inEvntOn((ttl > 0) & (inEvntOn > 0));%number of events when input TTL ports are in state 'ON'
        
        ttlPrtOn = unique(evntStr(inEvntOn));%different TTL ports (input ports only) in state 'ON'
        ttlIn(1:length(ttlPrtOn)) = struct('t', zeros(numel(inEvntOn), 2));%initialization
        origTTL(1:length(ttlPrtOn)) = struct('t', []);%initialization
        for ch = 1:length(ttlPrtOn) %run over different TTL ports
            z = 1;%counter of synchroimpulses
            for t = inEvntOn %run over inputs events when TTL set to On'
                if strcmp(evntStr{t}, ttlPrtOn{ch}) %right number of inputs port
                    ttlIn(ch).t(z, 1) = evntTStmp(t);%"on" (allStims(:, 1)) stimulus (mks from beginnig of day)
                    for n = inEvntOff(inEvntOff > t) %run over input events when TTL set to 'Off'
                        if strcmp(evntStr{n}(1:(end - 12)), ttlPrtOn{ch}(1:(end - 12))) %right number of inputs port
                            ttlIn(ch).t(z, 2) = evntTStmp(n);%"off"(allStims(:, 2)) stimulus (mks from beginnig of day)
                            z = z + 1;%counter of synchroimpulses
                            break;%out of (for n)
                        end
                    end
                end
            end
            ttlIn(ch).t(z:end, :) = [];%delete excess
            origTTL(ch).t = ttlIn(ch).t;%original timestamps of input TTLs
            
            for z = (numel(mStStRec) - 1):-2:2 %run over start-stop events
                jj = (ttlIn(ch).t(:, 1) >= mStStRec(z));%number of timestamps satisfying conditions
                ttlIn(ch).t(jj, :) = ttlIn(ch).t(jj, :) - (mStStRec(z) - mStStRec(z - 1)) + (512 * dsi);%stimulus moments from record begin
            end
            ttlIn(ch).t = ttlIn(ch).t - mStStRec(1);%adduction to zeros (first sample, mks from start-record-event)
            ttlIn(ch).t = ttlIn(ch).t / dsi;%convert stimulus moments to samples from record begin (from 0 because of arbitrary position of synchro-event in relation to sampling time)
            origTTL(ch).t = origTTL(ch).t - mStStRec(1);%adduction to zeros (first sample)
        end
    end
else
    disp('events-file does not exist')
    data = []; ttlIn = []; hd = []; spkTS = []; spkSM = [];
    return;%out of function
end
if isempty(strt)
    strt = 0;%begin read from recordation start (ms from record begin)
end
if (isempty(stp) || isequal(stp, 'e')) %read upto end of recordation
    stp = (diff(cscTStmp([1, end])) + (512 * dsi)) / 1e3;%read until end of recordation (ms)
end
rRecrd = [0, 0];%number of 512-samples-records to be read (initiation of array)

%= read data from ncs-files =%
data = zeros(1, rChNum);%LFP samples
if (nargout > 3)
    spkTS(1:rChNum, 1) = struct('tStamp', []);%spikes moment (mks)
end
if (nargout > 4)
    spkSM(1:rChNum, 1) = struct('shape', []);%spikes samples
end

n = 1;%number of channel in list channels to be read
csctstmpCopy = cscTStmp;%copy of all samples timestamps
for ch = rCh %run over all channels
    %%% lfp %%%
    fileToRead = ncsFiles(ch).f;%[pf, 'CSC', num2str(ch), '.ncs'];%pathname of file to be read
    if (ncsFiles(ch).r && exist(fileToRead, 'file')) %requested ncs-file exist and we can read file fully
        if isempty(hd) %no header
            cscHd = Nlx2MatCSC(fileToRead, [0 0 0 0 0], 1, 1, []);%header (lfp)
            adBitVolts = NlxParametr(cscHd, 'ADBitVolts');%multiplier to convert from samples to volts (lfp)
            dspDelay_mks = NlxParametr(cscHd, 'DspFilterDelay_탎');%DspFilterDelay_탎 (lfp)
            if ~isempty(dspDelay_mks)
                dspDelay_mks = dspDelay_mks * double(isequal(NlxParametr(cscHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (lfp)
            end
            inpInvert = double(strcmp(NlxParametr(cscHd, 'InputInverted'), 'True'));%input inverted
        else %we have a header
            adBitVolts = hd.adBitVolts(ch);%multiplier to convert from samples to volts (lfp)
            dspDelay_mks = hd.dspDelay_mks(ch);%DspFilterDelay_탎 (lfp)
            inpInvert = hd.inverted(ch);%input inverted
        end
        cscTStmp = csctstmpCopy;%Nlx2MatCSC(fileToRead, [1 0 0 0 0], 0, 1, []);%get timestamps of all samples
        %= get requested samples only =%
        rRecrd(1) = find(cscTStmp <= (cscTStmp(1) + (strt * 1e3)), 1, 'last');%number of first 512-samples-record to be read
        rRecrd(2) = find(cscTStmp <= (cscTStmp(1) + (stp * 1e3)), 1, 'last');%number of last 512-samples-record to be read
        cscTStmp = cscTStmp(rRecrd(1):rRecrd(2));%timestamps of requested 512-samples-record only
        smpl = Nlx2MatCSC(fileToRead, [0 0 0 0 1], 0, 2, rRecrd);%get requested samples
        %= end of get requested samples only =%

%         %= get all samples =%
%         smpl = Nlx2MatCSC(fileToRead, [0 0 0 0 1], 0, 1, []);%get all samples
%         %= end of get all samples =%

        if any(diff(cscTStmp, 2) > 5) %%(2)complex way (segment-mode or datasections lost; use timestamps)
            %disp('large error of timestamp')
            for z = (numel(mStStRec) - 1):-2:2 %run over start-stop events
                jj = (cscTStmp >= mStStRec(z));%number of timestamps satisfying conditions
                cscTStmp(jj) = cscTStmp(jj) - (mStStRec(z) - mStStRec(z - 1)) + (512 * dsi);%exclude interrecords periods
            end
            cscTStmp = cscTStmp - cscTStmp(1);%adduction the first sample to zeros

            data(floor(diff(cscTStmp([1, end])) / dsi), n) = 0;%memory allocation
            for t = 1:length(cscTStmp) %run over timestamps
                k = floor(cscTStmp(t) / dsi) + 1;
                data(k + (0:511), n) = smpl(:, t);%samples
            end
            data((k + 512):end, n) = [];%delete excess
        else %(1)simple way (concatenation)
            data(1:numel(smpl), n) = smpl(:);%samples
        end

        data(:, n) = data(:, n) * adBitVolts * 1e6;%microvoltes (lfp)
        if (inpInvert >= 1) %inverted signal
            data(:, n) = -1 * data(:, n);%back inverse
        end
    end
    
    %%% spikes %%%
    if (nargout > 3)
        fileToRead = [pf, 'SE', num2str(ch), '.nse'];%pathname of file to be read
        if exist(fileToRead, 'file') %requested file with spikes exist
            if isempty(hd) %no header
                spkHd = Nlx2MatSpike(fileToRead, [0 0 0 0 0], 1, 1, []);%header (spikes)
                adBitVoltsSpk = NlxParametr(spkHd, 'ADBitVolts');%multiplier to convert from samples to volts (spikes)
                dspDelay_mksSpk = NlxParametr(spkHd, 'DspFilterDelay_탎');%DspFilterDelay_탎 (spikes)
                dspDelay_mksSpk = dspDelay_mksSpk * double(isequal(NlxParametr(spkHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (spikes)
            else %we have a header
                adBitVoltsSpk = hd.adBitVoltsSpk(ch);%multiplier to convert from samples to volts (spikes)
                dspDelay_mksSpk = hd.dspDelay_mksSpk(ch);%DspFilterDelay_탎 (spikes)
                dspDelay_mks = hd.dspDelay_mks(ch);%DspDelayCompensation (scs)
            end
            spkTmStmp = Nlx2MatSpike(fileToRead, [1 0 0 0 0], 0, 1, []);%timestamps of spikes
%             spkTmStmp = spkTmStmp((spkTmStmp >= rStmps(1)) & (spkTmStmp <= rStmps(2)));%wanted spikes only
%                 %or:
%                 %spkTmStmp = Nlx2MatSpike(fileToRead, [1 0 0 0 0], 0, 4, rStmps);

            spkTS(n).tStamp = spkTmStmp - evntTStmp(1) - round((dspDelay_mksSpk + dspDelay_mks) / 2);%mks from record start
        end
        if (nargout > 4) %spikes time course requested
%             spkSM(n).shape = squeeze(Nlx2MatSpike(fileToRead, [0 0 0 0 1], 0, 4, rStmps));
            spkSM(n).shape = squeeze(Nlx2MatSpike(fileToRead, [0 0 0 0 1], 0, 1, []));
            spkSM(n).shape = spkSM(n).shape * adBitVoltsSpk * 1e6;%samples to microvolts
        end
    end
    n = n + 1;%number of channel in list channels to be read
end

%= cut nose and tail from data =%
%fileToRead = ncsFiles(rCh(1)).f;%[pf, 'CSC', num2str(ch), '.ncs'];%pathname of file to be read
cscTStmp = csctstmpCopy;%Nlx2MatCSC(fileToRead, [1 0 0 0 0], 0, 1, []);%get timestamps of all samples
rRecrd(1) = find(cscTStmp <= (cscTStmp(1) + (strt * 1e3)), 1, 'last');%number of first 512-samples-record to be read
rRecrd(2) = find(cscTStmp <= (cscTStmp(1) + (stp * 1e3)), 1, 'last');%number of last 512-samples-record to be read

cutSmplN = round(((cscTStmp(1) + (strt * 1e3)) - cscTStmp(rRecrd(1))) / dsi) - 1;%number of nose-samples to be cutted
if (cutSmplN > 0) %need to cut nose
    data(1:cutSmplN, :) = [];%cut nose
end
cutSmplN = 512 - round(((cscTStmp(1) + (stp * 1e3)) - cscTStmp(rRecrd(2))) / dsi) - 1;%number of tail-samples to be cutted
if (cutSmplN > 0) %need to cut tail
    data((end - cutSmplN):end, :) = [];%cut tail
end
%= end of cut nose and tail from data =%

    
%%% header compile (Axon Binary File (ABF) compatible)%%%
if (nargout > 2) %header requested
    hd.fFileSignature = 'Neuralynx';
    hd.nOperationMode = 3;%data were acquired in gap-free mode (continuous record)
    hd.lActualEpisodes = 1;%number of sweeps (for compatibility with abfload)

    hd.nADCNumChannels = 0;%number of channels
    dirCnt = dir(pf);%content of directory
    for t = 1:length(dirCnt)
        if ((~dirCnt(t).isdir) && (length(dirCnt(t).name) > 3))
            if isequal(dirCnt(t).name(end - 3:end), '.ncs')
                hd.nADCNumChannels = hd.nADCNumChannels + 1;%counter of channels
            end
        end
    end
    
    hd.adBitVolts = zeros(hd.nADCNumChannels, 1);%multiplier to convert from samples to volts (lfp)
    hd.dspDelay_mks = zeros(hd.nADCNumChannels, 1);%DspFilterDelay_탎 (lfp)
    hd.adBitVoltsSpk = zeros(hd.nADCNumChannels, 1);%multiplier to convert from samples to volts (spikes)
    hd.dspDelay_mksSpk = zeros(hd.nADCNumChannels, 1);%DspFilterDelay_탎 (spikes)
    hd.alignmentPt = zeros(hd.nADCNumChannels, 1);%spike samples back (from peak, including peak point)
    hd.inverted = zeros(hd.nADCNumChannels, 1);%input inverted
    hd.recChUnits = cell(hd.nADCNumChannels, 1);%mesurement units
    hd.recChNames = cell(hd.nADCNumChannels, 1);%name of channels
    hd.ch_si = zeros(hd.nADCNumChannels, 1);%sample interval (mks)

    %read headers
    for ch = 1:hd.nADCNumChannels
        %%% CSC(ncs)-files (lfp) %%%
        fileToRead = ncsFiles(ch).f;%[pf, 'CSC', num2str(ch), '.ncs'];%pathname of file to be read
        if exist(fileToRead, 'file') %requested file with lfp exist
            cscHd = Nlx2MatCSC(fileToRead, [0 0 0 0 0], 1, 1, []);%header (lfp)
            hd.ch_si(ch) = (1e6 / NlxParametr(cscHd, 'SamplingFrequency'));%sample interval (mks)
            hd.adBitVolts(ch) = NlxParametr(cscHd, 'ADBitVolts');%multiplier to convert from samples to volts (lfp)
            dspDelay_mks = NlxParametr(cscHd, 'DspFilterDelay_탎');%DspFilterDelay_탎 (lfp)
            if isempty(dspDelay_mks)
                dspDelay_mks = Inf;
            end
            hd.dspDelay_mks(ch) = dspDelay_mks;%DspFilterDelay_탎 (lfp)
            hd.dspDelay_mks(ch) = hd.dspDelay_mks(ch) * double(isequal(NlxParametr(cscHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (lfp)
            hd.recChUnits{ch} = '킮';%mesurement units
            hd.recChNames{ch} = NlxParametr(cscHd, 'AcqEntName');%name of channels
            hd.inverted(ch) = double(strcmp(NlxParametr(cscHd, 'InputInverted'), 'True'));%input inverted
        end

        %%% SE(nse)-files (spikes) %%%
        fileToRead = [pf, 'SE', num2str(ch), '.nse'];%pathname of file to be read
        if exist(fileToRead, 'file') %requested file with spikes exist
            spkHd = Nlx2MatSpike(fileToRead, [0 0 0 0 0], 1, 1, []);%header (spikes)
            hd.adBitVoltsSpk(ch) = NlxParametr(spkHd, 'ADBitVolts');%multiplier to convert from samples to volts (spikes)
            hd.dspDelay_mksSpk(ch) = NlxParametr(spkHd, 'DspFilterDelay_탎');%DspFilterDelay_탎 (spikes)
            hd.dspDelay_mksSpk(ch) = hd.dspDelay_mksSpk(ch) * double(isequal(NlxParametr(spkHd, 'DspDelayCompensation'), 'Disabled'));%DspDelayCompensation (spikes)
            hd.alignmentPt(ch) = NlxParametr(spkHd, 'AlignmentPt');%spike samples back (from peak, including peak point)
        end
    end
    hd.dataPtsPerChan = size(data, 1);%samples per channel
    hd.dataPts = hd.dataPtsPerChan * hd.nADCNumChannels;%total number of recorded samples
    hd.cscTStmp = cscTStmp;%timestamps of recorded samples (read largest file)
    hd.si = max(hd.ch_si);%sample interval (mks)
    hd.fADCSampleInterval = hd.si;%sample interval (mks)

    %fill hd.recTime ('ttl' and 'evntStr' automatically exist)
    %if (~exist('ttl', 'var') || ~exist('evntStr', 'var'))
    %   fileToRead = [pf, 'Events.nev'];%pathname of file to be read
    %   [ttl, evntStr] = Nlx2MatEV(fileToRead, [0 0 1 0 1], 0, 1, []);%read event-file
    %end
    hd.inTTL_timestamps = origTTL;%original timestamps of input TTLs (adducted to zeros - first sample is zero)
    hd.TTLs = ttl;%ttl events
    hd.EventStrings = evntStr;%text of events
    if exist('cscHd', 'var')
        hd.cscHd = cscHd;%original header of neuralynx-file
    end
    if exist('spkHd', 'var')
        hd.spkHd = spkHd;%original header of neuralynx-file
    end

    fid = fopen([pf, 'CheetahLogFile.txt']);%open cheetah log-file
    logy = fread(fid, 'char');%read log-file
    logy = char(logy');
    fclose(fid);

    jj = strfind(logy, '*-');
    tmS = cell(numel(jj), 3);
    if (nlxVer == 0) %initial vertion of cheetah logs
        for t = 1:numel(jj)
            z = strfind(logy(jj(t):(jj(t) + 50)), ' - ');
            bufStr = textscan(logy((jj(t) + 2):(jj(t) + z(2) - 2)), '%s%s%s');
            tmS(t, 1) = bufStr{1}(1);%time of current event (hh:mm:ss.ms)
            tmS{t, 2} = str2double(bufStr{3}{1});%timestamp of current event (microseconds from Neuralynx on)
            tmS{t, 3} = logy((jj(t) + z(2) + 2):(jj(t) + z(2) + 2 + 33));%current event string
        end
        ymdhms = datevec(tmS{1, 1}, 'HH:MM:SS.FFF');%date as numels
    else %next vertion of cheetah logs
        for t = 1:numel(jj)
            z = strfind(logy(jj(t):(jj(t) + 50)), ' - ');
            tmS{t, 1} = logy((jj(t) + 2):(jj(t) + z(1) - 2));%time of current event (hh:mm:ss.ms)
            tmS{t, 2} = str2double(logy((jj(t) + 2 + z(1)):(jj(t) + z(2) - 2)));%timestamp of current event (microseconds from Neuralynx on)
            tmS{t, 3} = logy((jj(t) + z(2) + 2):(jj(t) + z(2) + 2 + 33));%current event string
        end
        ymdhms = datevec(tmS{1, 1}, 'yyyy/mm/dd HH:MM:SS');%date as numels
    end
    nlxBeginS = ymdhms(4) * 3600 + ymdhms(5) * 60 + ymdhms(6);%time of Neuralynx turned ON (seconds from day begin)
    
    n = 0;%number of first event 'Starting Recording'
    for t = 1:length(evntStr)
        if isequal(evntStr{t}, 'Starting Recording')
            n = t;%number of first event 'Starting Recording'
            break;%out of "for t"
        end
    end
    [~, t] = min(abs(vertcat(tmS{:, 2}) - evntTStmp(n)));%nearest timestamp
    jj = false;
    z = -1;
    while (~jj)
        z = z + 1;
        jj = strcmp(tmS{t + z, 3}, 'AcquisitionControl::StartRecording');%go forward to find
        if ~jj
            jj = strcmp(tmS{t - z, 3}, 'AcquisitionControl::StartRecording');%go backward to find
            if jj
                z = -z;
            end
        end
    end
    t = t + z;%number of row with pure time of start
    if (nlxVer == 0) %initial vertion of cheetah logs
        rcTm = datenum(['2011:01:01 ', tmS{t, 1}], 'yyyy:mm:dd HH:MM:SS.FFF') - datenum('2011:01:01 00:00:00.000', 'yyyy:mm:dd HH:MM:SS.FFF');%days after midnight
        hd.recTime(1) = rcTm * (24 * 60 * 60);%record stop time (seconds after midnight)
        if ((hd.recTime(1) - nlxBeginS) <= 0) %new day began
            hd.recTime(1) = hd.recTime(1) + (24 * 60 * 60);%record stop time (seconds after midnight)
        end
        %hd.recTime(1) = (tmS{t, 2}  / 1e6) + nlxBeginMKS;%record start time (seconds from day begin (seconds after midnight))
    else
        rcTm = datenum(tmS{t, 1}, 'yyyy/mm/dd HH:MM:SS') - datenum([tmS{1, 1}(1:11), ' 00:00:00'], 'yyyy/mm/dd HH:MM:SS');%days after midnight
        hd.recTime(1) = rcTm * (24 * 60 * 60);%record stop time (seconds after midnight)
    end
    
    n = 0;%number of last event 'Stopping Recording'
    for t = length(evntStr):-1:1
        if isequal(evntStr{t}, 'Stopping Recording')
            n = t;%number of last event 'Stopping Recording'
            break;%out of "for t"
        end
    end
    [~, t] = min(abs(vertcat(tmS{:, 2}) - evntTStmp(n)));%nearest timestamp
    jj = false;
    z = -1;
    while (~jj)
        z = z + 1;
        jj = strcmp(tmS{t + z, 3}, 'AcquisitionControl::StopRecording(');
        if ~jj
            jj = strcmp(tmS{t - z, 3}, 'AcquisitionControl::StopRecording(');
            if jj
                z = -z;
            end
        end
    end
    t = t + z;%number of row with pure time of stop
    if (nlxVer == 0) %initial vertion of cheetah logs
        rcTm = datenum(['2011:01:01 ', tmS{t, 1}], 'yyyy:mm:dd HH:MM:SS.FFF') - datenum('2011:01:01 00:00:00.000', 'yyyy:mm:dd HH:MM:SS.FFF');%days after midnight
        hd.recTime(2) = rcTm * (24 * 60 * 60);%record stop time (seconds after midnight)
        if ((hd.recTime(2) - nlxBeginS) <= 0) %new day began
            hd.recTime(2) = hd.recTime(2) + (24 * 60 * 60);%record stop time (seconds after midnight)
        end
    else
        rcTm = datenum(tmS{t, 1}, 'yyyy/mm/dd HH:MM:SS') - datenum([tmS{1, 1}(1:11), ' 00:00:00'], 'yyyy/mm/dd HH:MM:SS');%days after midnight
        hd.recTime(2) = rcTm * (24 * 60 * 60);%record stop time (seconds after midnight)
    end
    
    %(hd.sweepStartInPts * hd.fADCSampleInterval)   the start times of sweeps in sample points (from beginning of recording)
    %hd.sweepStartInPts = ?allStims(:, 1)?;%the start times of sweeps in sample points (from beginning of recording)
end

function mStStRec = FindStStStemp(cscTS, evntStr, ststTStmp, si)
%find moments of "Starting Recording" and "Stopping Recording"
%
%INPUTS
%cscTS - timestamps of recorded samples (mks, largest file)
%evntStr - strings with events description
%ststTStmp - timestampes of start-stop and ttl events
%si - sample interval (mks)
%
%OUTPUTS
%mStStRec - numbers of timestamp of start and stop recordings
%

%methode 1 (find by events timestamp)
mStStRec = zeros(length(evntStr), 1);%preallocation of memory
z = 1;
for t = 1:length(evntStr)
    if strcmp(evntStr{t}, 'Starting Recording')
        mStStRec(z) = ststTStmp(t);%start timestamp (number)
        z = z + 1;
    end
    if strcmp(evntStr{t}, 'Stopping Recording')
        mStStRec(z) = ststTStmp(t);%stop timestamp (number)
        z = z + 1;
    end
end
mStStRec(z:end, :) = [];%delete excess
for t = 1:length(mStStRec)
    [~, z] = min(abs(cscTS - mStStRec(t)));%number of nearest timestemp of samples-block
    ststrec = cscTS(z);%nearest timestemp of samples-block
    while (abs(ststrec - mStStRec(t)) > (512 * si)) %lost records detected on edges
        if ((t / 2) == round(t / 2)) %even timestamp (stop-event)
            ststrec = ststrec + (512 * si);%compensation of lost record
        else %odd timestamp (start-event)
            ststrec = ststrec - (512 * si);%compensation of lost record
        end
    end
    mStStRec(t) = ststrec;%replace number of timestamp by corresponding timestamp
end

%methode 2 (find by lfp timestamp)
% difTS = diff(evntTS);%difference
% tmp = find(difTS > ((512 * si) + 10));%number of timestamp with "Stop" events
% mStStRec = zeros((2 * numel(tmp)) + 2, 1);%preallocation of memory
% mStStRec(2:2:(end - 2)) = tmp;%stop timestamps
% mStStRec(3:2:(end - 1)) = tmp + 1;%start timestamps
% mStStRec(1) = 1;%first timestamp corresponds to first start
% mStStRec(end) = length(evntTS);%last timestamp corresponds to last stop

function nlxPrm = NlxParametr(headCell, fieldNm)
%get value of specified parameter
%
%INPUTS
%cscHd - cell array wiht Neuralynx parameter
%fieldNm - name of requested parameter (single name)
%
%OUTPUTS
%nlxPrm - value of parameter

nlxPrm = [];%initialization
for n = 1:length(headCell) %run over cells with parameters
    if ~isempty(strfind(headCell{n}, fieldNm)) %string contains requested name
        t = find(headCell{n} == ' ', 1, 'first');%find delimiter
        strVal = headCell{n}((t + 1):end);%string with value of parameter
        if (any(double(strVal) < 46) || any(double(strVal) > 57)) %value is word
            nlxPrm = strVal;
        else %value is numeric
            nlxPrm = str2double(strVal);%numeric value of parameter
        end
        break;%out of (for n)
    end
end
