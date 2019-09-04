function [data, synchro, hd, spkTS, spkSM] = ZavLoadData(flNm, inHd, rCh, rSw, strt, stp, nlxVer)
%[data, hd, synchro, spkTS] = ZavLoadData(flNm, inHd, rCh, rSw, strt, stp)
%unified interface for load files (abf-files, Neuralynx, edf+)
%
%INPUTS
%flNm - pathname of file
%inHd - input head
%rCh - numbers of channels to be read
%rSw - numbers of sweeps to be read
%strt - start time (ms from record begin)
%stp - stop time (ms from record begin)
%nlxVer - version ov NLX data (0 - initial, 1 - new, ...)
%
%OUTPUTS
%data - signal samples
%synchro - stimulus moments (samples from sweep beginning)
%hd - file header (information about record)
%spkTS - spikes appearence moments (mks)
%spkSM - spikes samples (Neuralynx)
%

t = find(flNm(1:end) == '.', 1, 'last');
if (~isempty(t) && ((length(flNm) - t) < 7))
    ext = flNm(t:end);%file extention exist
else
    ext = '';%no file extention
end

if (nargout > 1) %synchro requested
    synchro = struct('t', []);%initial empty
end
if (nargout > 3) %spikes timestampes requested
    spkTS = [];%initial empty
end
if (nargout > 4) %spikes samples requested
    spkSM = [];%initial empty
end

switch ext
    case '.abf' %axon binary file
        if isempty(rCh)
            rCh = 'a';%read all channels
        end
        if isempty(rSw)
            rSw = 'a';%read all sweeps
        end
        if isempty(strt)
            strt = 0;%begin read from recordation start
        end
        strt = strt / 1e3;%convert to seconds
        if isempty(stp)
            stp = 'e';%read until end of recordation
        else
            stp = stp / 1e3;%convert to seconds
        end
        
        if (~isequal(rCh, 'a') && isnumeric(rCh))
            if isempty(inHd)
                [~, ~, inHd] = abfload(flNm, 'channels', 'a', 'sweeps', 1, 'start', 0, 'stop', 0.1);%get header
            end
            abfCh = cell(1, length(rCh));%names of requested channels
            for t = 1:length(rCh)
                abfCh{t} = inHd.recChNames{rCh(t)};%names of requested channels
            end
        else %rCh containes channel names (must be cell or string) or string 'a'
            abfCh = rCh;%names of rquested channels
        end
        [data, ~, hd] = abfload(flNm, 'channels', abfCh, 'sweeps', rSw, 'start', strt, 'stop', stp);
        if (hd.nOperationMode == 3) %gap-free mode
            hd.lActualEpisodes = 1;%number of sweeps
        end
    case {'.ncs', '.nse', '.nev', ''} %Neuralynx
        if isempty(nlxVer)
            nlxVer = 0;%old version of Neuralynx Cheetah
        end
        switch nargout
            case 5 %all variables requested
                [data, synchro, hd, spkTS, spkSM] = ZavNrlynx(flNm, inHd, rCh, nlxVer, strt, stp);%load 
            case 4 %first 4 variables requested
                [data, synchro, hd, spkTS] = ZavNrlynx(flNm, inHd, rCh, nlxVer, strt, stp);%load 
            case 3 %first 3 variables requested
                [data, synchro, hd] = ZavNrlynx(flNm, inHd, rCh, nlxVer, strt, stp);%load 
            case 2 %first 2 variables requested
                [data, synchro] = ZavNrlynx(flNm, inHd, rCh, nlxVer, strt, stp);%load 
            case 1 %samples requested only
                data = ZavNrlynx(flNm, inHd, rCh, nlxVer, strt, stp);%load 
        end
    case '.daq' %daq-files;
        switch nargout
            case {3, 4, 5}
                [data, hd] = ZavDaqRead(flNm, inHd, rCh, rSw, strt, stp);%read specified channels
            case {1, 2}
                data = ZavDaqRead(flNm, inHd, rCh, rSw, strt, stp);%read specified channels
        end
    case '.edf' %eeg in "european data format"
        [data, hd] = ZavManEEGload(flNm, rCh);%load EEG in "european data format"
    otherwise
        disp('error file type')
end
