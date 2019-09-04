function [data, hd] = ZavDaqRead(flNm, inHd, rCh, rSw, strt, stp)
%[data, hd] = ZavDaqRead(flNm, inHd, rCh, strt, stp)
%read specified channels from daq-files
%
%INPUTS
%flNm - pathname of file
%inHd - input head
%rCh - numbers of channels to be read
%rSw - numbers of sweeps to be read
%strt - start time (for 'abf' and 'daq' in seconds)
%stp - stop time (for 'abf' and 'daq' in seconds)
%
%OUTPUTS
%data - recorded samples (mkV) %synchro - moments of synchroevents (samples from record start)
%hd - information of header

if ((nargout > 1) || isempty(inHd))
    hd = daqread(flNm, 'info');%recordation attributes
elseif ~isempty(inHd)
    hd = inHd;%header
end
if (nargout > 1) %header requested
    hd.fFileSignature = 'DAQ';%type of file
    hd.lActualEpisodes = hd.ObjInfo.TriggersExecuted - 1;%number of sweeps
    hd.si = 1e6 / hd.ObjInfo.SampleRate;%sample interval (mks)
    hd.nADCNumChannels = length(hd.ObjInfo.Channel);%number of channels
    hd.nOperationMode = 1 + (2 * double(hd.lActualEpisodes < 1));%type of record
    dayTm = datenum([hd.ObjInfo.EventLog(1).Data.AbsTime(1:3), 0, 0, 0]);%time of day begin
    hd.recTime(1) = (datenum(hd.ObjInfo.EventLog(1).Data.AbsTime) - dayTm) * 3600 * 24;%start (seconds after day begin)
    hd.recTime(2) = (datenum(hd.ObjInfo.EventLog(end).Data.AbsTime) - dayTm) * 3600 * 24;%stop (seconds after day begin)
    hd.sweepLengthInPts = hd.ObjInfo.SamplesPerTrigger;
end
    
%     synchro = struct('t', zeros(hd.lActualEpisodes, 2));%preallocation of memory for synchroevents
%     n = 1;%counter of trigger
%     for t = 1:length(hd.ObjInfo.EventLog)
%         if isequal(hd.ObjInfo.EventLog(t).Type, 'Trigger')
%             synchro(1).t(n, 1:2) = hd.ObjInfo.EventLog(t).Data.RelSample;%moment of synchroevent (samples from record start)
%             n = n + 1;%counter of trigger
%         end
%         if (n > hd.lActualEpisodes)
%             break;%out of (for t)
%         end
%     end

if isequal(rCh, 'a') %all channels requested
    rCh = 1:hd.nADCNumChannels;%read all channels
end
if isequal(rSw, 'a')
    rSw = 1:hd.lActualEpisodes;%read all sweeps
end
if ~isempty(strt)%start time (for 'abf' and 'daq' in seconds)
    strt = round(strt * hd.ObjInfo.SampleRate);%start sample
    if (strt <= 0)
        strt = 1;%first sample
    end
else
    strt = 1;%first sample
end
if ~isempty(stp)%stop time (for 'abf' and 'daq' in seconds)
    if isequal(stp, 'e')
        stp = hd.ObjInfo.SamplesAcquired;%last sample
    else
        stp = round(stp * hd.ObjInfo.SampleRate);%stop sample
    end
else
    stp = hd.ObjInfo.SamplesAcquired;%last sample
end

if (isfinite(hd.ObjInfo.SamplesPerTrigger) && (numel(rSw) > 1)) %sweeps mode
    if (~isempty(diff(rSw)) && all(diff(rSw) == 1)) %number in a row
        data = daqread(flNm, 'Channels', rCh, 'Triggers', rSw([1, end]));%read specified channels
        data(isnan(data(:, 1)), :) = [];%exclude NaNs
        data = reshape(data, hd.ObjInfo.SamplesPerTrigger, numel(rSw), numel(rCh));%reshape
        data = permute(data, [1, 3, 2]);%reorder dimentions
    else %arbitrary order of sweep
        data = zeros(hd.ObjInfo.SamplesPerTrigger, numel(rCh), numel(rSw));%memory preallocation
        for t = 1:numel(rSw) %run over sweeps
            dataOrig = daqread(flNm, 'Channels', rCh, 'Triggers', rSw(t));%read specified channels
            for n = 1:numel(rCh)
                data(:, n, t) = dataOrig(:, n);
            end
        end
    end
else %continuous mode (gap free, single sweep)
    data = daqread(flNm, 'Channels', rCh, 'Samples', [strt, stp]);
end

%convert to microvolts
for t = 1:numel(rCh)
    ch = rCh(t);
    offSet = hd.ObjInfo.Channel(ch).NativeOffset;
    scal = hd.ObjInfo.Channel(ch).NativeScaling;
    data(:, t) = (data(:, t) - offSet) / scal;
    scal = 1;
%     switch hd.ObjInfo.Channel(ch).Units
%         case 'V'
%             scal = 1e6;
%         case 'mV'
            scal = 1e3;
%     end
    data(:, t, :) = data(:, t, :) * scal;%convert to mkV
end
