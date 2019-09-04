function varargout = ExpressAnalysis(varargin)
% EXPRESSANALYSIS MATLAB code for ExpressAnalysis.fig
%      EXPRESSANALYSIS, by itself, creates a new EXPRESSANALYSIS or raises the existing
%      singleton*.
%
%      H = EXPRESSANALYSIS returns the handle to a new EXPRESSANALYSIS or the handle to
%      the existing singleton*.
%
%      EXPRESSANALYSIS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EXPRESSANALYSIS.M with the given input arguments.
%
%      EXPRESSANALYSIS('Property','Value',...) creates a new EXPRESSANALYSIS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ExpressAnalysis_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ExpressAnalysis_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ExpressAnalysis

% Last Modified by GUIDE v2.5 04-Jun-2019 15:32:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ExpressAnalysis_OpeningFcn, ...
                   'gui_OutputFcn',  @ExpressAnalysis_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ExpressAnalysis is made visible.
function ExpressAnalysis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ExpressAnalysis (see VARARGIN)

% Choose default command line output for ExpressAnalysis
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes ExpressAnalysis wait for user response (see UIRESUME)
% uiwait(handles.ExpressAnalysis);


% --- Outputs from this function are returned to the command line.
function varargout = ExpressAnalysis_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

expAn = handles.ExpressAnalysis;%handle of main form (figure)
setappdata(expAn, 'pth', '');%path to file
setappdata(expAn, 'flNm', '');%name of file
setappdata(expAn, 'scanDrct', '');%merge directory
set(handles.ExportData, 'UserData', []);%null time (injection)



% --------------------------------------------------------------------
function ScanDir_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ScanDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%scan directory, file conversion
expAn = handles.ExpressAnalysis;
textMsgs = handles.TextMsgs;
scanDrct = getappdata(expAn, 'scanDrct');%merge directory
if isempty(scanDrct)
    scanDrct = getappdata(expAn, 'pth');%path
end
scanDrct = uigetdir(scanDrct, 'Directory with file to be converted');%get directory
setappdata(expAn, 'scanDrct', scanDrct);
scanDrct = [scanDrct, '\'];
dirCont = dir(scanDrct);%content of directory

set(textMsgs, 'String', {});
clear fBase;
fBase(1:1000) = struct('srcFile', '', 'cscN', [], 'chnlGrp', [], 'type', '', 'matFile', '');
w = 1;%files counter
ii = zeros(1000, 1);%file order
for t = 3:length(dirCont) %run over directory content
    if (dirCont(t).isdir) %NLX
        if exist([scanDrct, dirCont(t).name, '\Events.nev'], 'file')
            fBase(w).srcFile = [scanDrct, dirCont(t).name, '\'];%Nlx
            fBase(w).type = 'NLX';%type of file
            %ttl = Nlx2MatEV([fBase(w).srcFile, 'Events.nev'], [0 0 1 0 0], 0, 1, []);
            subDirCnt = dir(fBase(w).srcFile);
            fBase(w).cscN = 0;%number of LFP channels
            for z = 3:length(subDirCnt) %run over objects of current folder
                if (~subDirCnt(z).isdir && (length(subDirCnt(z).name) > 3)) %file with quite long name
                    if (isequal(subDirCnt(z).name((end - 3):end), '.ncs') && ~isempty(strfind(subDirCnt(z).name, 'CSC'))) %ncs channel (LFP trace)
                        fBase(w).cscN = fBase(w).cscN + 1;%number of LFP channels
                    end
                end
            end
            
            %find config file (read log)
            fid = fopen([fBase(w).srcFile, 'CheetahLogFile.txt'], 'r');
            clgfl = char(fread(fid)');
            fclose(fid);
            cfgIn = strfind(clgfl, 'NexusConfigs\');%starting indices of config file extensions
            if ~isempty(cfgIn) %subdirectory exist
                bufStr = clgfl(cfgIn(end) + (13:60));%part of configfile pathname
                z = strfind(bufStr, '.cfg');
                cfgfl = bufStr(1:(z(1) + 3));
                fBase(w).chnlGrp = cell(0);
                if isequal(cfgfl, 'A1x16_Brd1_double_4piezo_7emg_termo.cfg')
                    fBase(w).chnlGrp(1) = {1:16}; fBase(w).chnlGrp(2) = {17:32}; fBase(w).chnlGrp(3) = {33:44};
                elseif isequal(cfgfl, 'A1x16_Brd1_double_4piezo_4emg_4DCdiff.cfg')
                    fBase(w).chnlGrp(1) = {1:16}; fBase(w).chnlGrp(2) = {17:32}; fBase(w).chnlGrp(3) = {33:44};
                elseif isequal(cfgfl, 'A1x32_Brd1_double_4piezo_4emg_4DCdiff.cfg')
                    fBase(w).chnlGrp(1) = {1:32}; fBase(w).chnlGrp(2) = {33:64}; fBase(w).chnlGrp(3) = {65:76};
                elseif isequal(cfgfl, 'A1x32_Brd1_double_4piezo_4emg_2DCdiff.cfg')
                    fBase(w).chnlGrp(1) = {1:32}; fBase(w).chnlGrp(2) = {33:64}; fBase(w).chnlGrp(3) = {65:74};
                elseif isequal(cfgfl, 'A1x32_Brd2_4piezo_4emg_4DCdiff.cfg')
                    fBase(w).chnlGrp(1) = {1:32}; fBase(w).chnlGrp(2) = {33:44};
                elseif isequal(cfgfl, 'A1x16_Brd1_double_piezo_ecg_ZAV_diff.cfg')
                    fBase(w).chnlGrp(1) = {1:16}; fBase(w).chnlGrp(2) = {17:32}; fBase(w).chnlGrp(3) = {33:36};
                elseif isequal(cfgfl, 'A1x16_Brd1_double_4piezo.cfg')
                    fBase(w).chnlGrp(1) = {1:16}; fBase(w).chnlGrp(2) = {17:32}; fBase(w).chnlGrp(3) = {33:36};
                elseif isequal(cfgfl, 'A1x16_Brd2_double_4piezo.cfg')
                    fBase(w).chnlGrp(1) = {1:16}; fBase(w).chnlGrp(2) = {17:32}; fBase(w).chnlGrp(3) = {33:36};
                end
            end
            
            ii(w) = t;%files order
            w = w + 1;%next file
        end
    elseif ((length(dirCont(t).name) > 3) && isequal(dirCont(t).name((end - 3):end), '.daq')) %daq
        fBase(w).srcFile = [scanDrct, dirCont(t).name];%daq
        fBase(w).type = 'DAQ';%type of file
        fBase(w).cscN = NaN;
        ii(w) = t;%files order
        w = w + 1;%next file
    elseif ((length(dirCont(t).name) > 3) && isequal(dirCont(t).name((end - 3):end), '.abf')) %abf
        fBase(w).srcFile = [scanDrct, dirCont(t).name];%abf
        fBase(w).type = 'ABF';%type of file
        fBase(w).cscN = NaN;
        ii(w) = t;%files order
        w = w + 1;%next file
    end
end
for t = 3:length(dirCont) %run over directory content
    if ((length(dirCont(t).name) > 3) && isequal(dirCont(t).name((end - 3):end), '.mat')) %abf
        bufStr = [scanDrct, dirCont(t).name(1:(end - 4))];%base name
        for z = 1:length(fBase)
            if strfind(fBase(z).srcFile, bufStr) %source file exist
                fBase(z).matFile = [bufStr, '.mat'];%add field
                bufStr = '';%file accepted
                break;%out of (for z)
            end
        end
        if ~isempty(bufStr) %original mat-file
            if ((length(bufStr) < 5) || isempty(strfind(bufStr((end - 4):end), '_brst')))
                fBase(w).matFile = [bufStr, '.mat'];
                ii(w) = t;%files order
                w = w + 1;%next file
            end
        end
    end
end
fBase(w:end) = [];%delete empty
ii(w:end) = [];%delete excess
[~, ii] = sort(ii);
fBase = fBase(ii);

msgs = cell(length(fBase), 1);
rList = cell(length(fBase), 1);%list of recordations to be treated
z = length(num2str(length(fBase)));
for t = 1:length(fBase)
    bufStr = num2str(t);
    bufStr = [repmat(' ', 1, z), bufStr];
    bufStr = bufStr((end - z + 1):end);
    if ~isempty(fBase(t).srcFile)
        if isequal(fBase(t).type, 'NLX')
            k = find(fBase(t).srcFile(1:(end - 2)) == '\', 1, 'last') + 1;
        else
            k = find(fBase(t).srcFile == '\', 1, 'last') + 1;
        end
        msgs{t} = [bufStr, ': ', fBase(t).srcFile(k:end)];
        if ~isempty(fBase(t).matFile)
            msgs{t} = [msgs{t}, ' (+mat)'];
        end
    else
        k = find(fBase(t).matFile == '\', 1, 'last') + 1;
        msgs{t} = [bufStr, ': ', fBase(t).matFile(k:end)];
    end
    rList{t, 1} = msgs{t}; rList{t, 2} = true;
end
set(textMsgs, 'String', msgs);
set(handles.recList, 'Data', rList);

setappdata(expAn, 'fBase', fBase);%set list of recordations
set(handles.ExportData, 'Enable', 'on');%figures plotting allowed

%erase data of single recordation
setappdata(expAn, 'hd', [])
setappdata(expAn, 'lfp', [])
setappdata(expAn, 'lfpVar', [])
setappdata(expAn, 'spks', [])
setappdata(expAn, 'zavp', [])



% --- Executes on button press in ConvertList.
function ConvertList_Callback(hObject, eventdata, handles)
% hObject    handle to ConvertList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%analysis of data
expAn = handles.ExpressAnalysis;%main form handle
textMsgs = handles.TextMsgs;%messeges
fBase = getappdata(expAn, 'fBase');%list of files
rList = get(handles.recList, 'Data');%list of wanted files
rList = vertcat(rList{:, 2});%flags of accepted files
for rf = 1:length(fBase) %run over records
    if rList(rf) %treat wanted files only
        if isempty(fBase(rf).matFile) %non converted data
            if isequal(fBase(rf).type, 'NLX') %NLX
                nlxVer = 0;%old format of Cheetah log-file
                dirCnt = dir(fBase(rf).srcFile);%directory content
                for z = 3:length(dirCnt) %run over files in NLX directory
                    if ((length(dirCnt(z).name) > 3) && isequal(dirCnt(z).name(end - 3:end), '.nde'))
                        nlxVer = 1;%new format of Cheetah log-file
                        break;%out of (for z)
                    end
                end
                ttl = Nlx2MatEV([fBase(rf).srcFile, 'Events.nev'], [0 0 1 0 0], 0, 1, []);%stimulus pattern
                try
                    ZavSoursAnalys(fBase(rf).srcFile, 1, any(ttl > 0), [], 1:fBase(rf).cscN, fBase(rf).chnlGrp, nlxVer);%convertation
                    fBase(rf).matFile = [fBase(rf).srcFile(1:(end - 1)), '.mat'];%corresponding mat-file
                catch errMess
                    msgs = get(textMsgs, 'String');
                    msgs{end + 1} = ['error on ', fBase(rf).srcFile, ' :=:', errMess.message];
                    set(textMsgs, 'String', msgs);
                end
            else %DAQ or ABF
                try
                    isStim = false;%evoked or spontaneous activity
                    ZavSoursAnalys(fBase(rf).srcFile, 1, isStim, [], [], [], 0);%convertation
                    fBase(rf).matFile = [fBase(rf).srcFile(1:(end - 4)), '.mat'];%corresponding mat-file
                catch errMess
                    msgs = get(textMsgs, 'String');
                    msgs{end + 1} = ['error on ', fBase(rf).srcFile, ' :=:', errMess.message];
                    set(textMsgs, 'String', msgs);
                end
            end
        end
        msgs = get(textMsgs, 'String');
        msgs{end + 1} = [num2str(rf), ' of ', num2str(length(fBase)), ' done (', num2str(fBase(rf).cscN), ' csc)'];
        set(textMsgs, 'String', msgs);

        %treatAllRecs = isequal(get(handles.TreatAllRecords, 'State'), 'on');%treat all records
        if (isequal(get(handles.TreatAllRecords, 'State'), 'on') && ~isempty(fBase(rf).matFile)) %treat all records
            %send parameters of current file
            k = find(fBase(rf).matFile == '\', 1, 'last');
            pth = fBase(rf).matFile(1:k);%directory
            flNm = fBase(rf).matFile((k + 1):end);%filename
            setappdata(expAn, 'pth', pth);%path
            setappdata(expAn, 'flNm', flNm);%file
            if isempty(get(handles.PrincChan, 'String')) %request principal channel
                set(handles.PrincChan, 'String', '1')
            end
            varInFl = who(matfile([pth, flNm]));%variables in requested file
            load([pth, flNm])
            if ismember('zp', varInFl)
                zavp = zp;
            end
            setappdata(expAn, 'hd', hd)
            setappdata(expAn, 'lfp', double(lfp))
            setappdata(expAn, 'lfpVar', lfpVar)
            setappdata(expAn, 'spks', spks)
            setappdata(expAn, 'zavp', zavp)
            if exist('sepOnsetPeak', 'var') %variable exist
                setappdata(expAn, 'sepOnsetPeak', sepOnsetPeak);%SEP parameters
            elseif isappdata(expAn, 'sepOnsetPeak') %previouse data 
                rmappdata(expAn, 'sepOnsetPeak');%remove previouse data
            end
            AnalyseSinglEvoked(handles);%calculcate SEP parameters
        end
    end
end
setappdata(expAn, 'fBase', fBase);%set list of recordations



% --------------------------------------------------------------------
function OpenFile_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to OpenFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

expAn = handles.ExpressAnalysis;%main form handle
pth = getappdata(expAn, 'pth');%directory
flNm = getappdata(expAn, 'flNm');%filename
if isempty(pth) %no file choosed
    flNm = '';
    pth = getappdata(expAn, 'scanDrct');
end

[flNm, pth] = uigetfile({'.mat'}, 'Select file', [pth, flNm]);%open dialog
if (isequal(flNm, 0) && isequal(pth, 0)) %no file choosed
    return
end
setappdata(expAn, 'pth', pth);%path
setappdata(expAn, 'flNm', flNm);%file

if isempty(get(handles.PrincChan, 'String')) %request principal channel
    errMsgH = errordlg('Enter single channel as principal and click ''Analyse''', 'No channel specified', 'modal');
    return;%do not analysis
end
load([pth, flNm])
setappdata(expAn, 'hd', hd)
setappdata(expAn, 'lfp', double(lfp))
setappdata(expAn, 'lfpVar', lfpVar)
setappdata(expAn, 'spks', spks)
setappdata(expAn, 'zavp', zavp)
if exist('sepOnsetPeak', 'var') %variable exist
    setappdata(expAn, 'sepOnsetPeak', sepOnsetPeak);%SEP parameters
elseif isappdata(expAn, 'sepOnsetPeak') %previouse data 
    rmappdata(expAn, 'sepOnsetPeak');%remove previouse data
end
% setappdata(expAn, 'fBase', []);%erase list of recordations
AnalyseSinglEvoked(handles);%perform analysis



function AnalyseSinglEvoked(handles)
%analysis of experimental data
% handles    structure with handles and user data (see GUIDATA)

expAn = handles.ExpressAnalysis;%main form handle
textMsgs = handles.TextMsgs;%messeges handle
hd = getappdata(expAn, 'hd');%header
if isempty(hd) %no recordation opened
    %set(textMsgs, 'String', {'Choose directory', 'or', 'Open mat-file'});
    return
end


% %load common-zero
% load('2017-12-08_wLFPcomm')
% wLFPcomm = wCmn;%mean(wLFPcomm(:, :, 1:end), 3);


ax1 = handles.Axes1;%handle of axes
flNm = getappdata(expAn, 'flNm');%filename
set(textMsgs, 'String', {});
if isequal(flNm((end - 3):end), '.mat')
    zavp = getappdata(expAn, 'zavp');%load([pth, flNm], 'zavp')
    
    %ON-response of evoked signals
    if ~isempty(zavp.realStim(1).r) %evoked signals (on-response)
        pth = getappdata(expAn, 'pth');%directory
        %hd = getappdata(expAn, 'hd');%header
        lfp = getappdata(expAn, 'lfp');%LFP
        zavp = getappdata(expAn, 'zavp');%metadata
        treatAllRecs = isequal(get(handles.TreatAllRecords, 'State'), 'on');%treat all records
        treatAllCh = isequal(get(handles.TreatAllChannels, 'State'), 'on');%treat all channel
        if isappdata(expAn, 'sepOnsetPeak') %variable exist
            sepOnsetPeak = getappdata(expAn, 'sepOnsetPeak');%SEP parameters
        end
        set(expAn, 'Name', ['Express analysis - ', pth, flNm]);%full name
        segms = GetSegments(hd, zavp, 1);%segments of evoked activity
        setappdata(expAn, 'segms', segms);%segments
        rCh = str2double(get(handles.PrincChan, 'String'));%zavp.chL4(1);%principal channel
        setappdata(expAn, 'rCh', rCh);%principal channel

        ii = [];%principal channels
        if exist('sepOnsetPeak', 'var') %SEP parameters was calculated earlier
            [sepOnsetPeak, ii] = CorrectSepOP(expAn, sepOnsetPeak, hd, zavp);%correct sepOnsetPeak
            sopExist = true;%flag of existing variable
        end
        if (~exist('sepOnsetPeak', 'var') || ~ismember(rCh, ii) || treatAllCh) %no variable or new principal channel
            sopExist = false;%flag of existing variable
            if ~exist('sepOnsetPeak', 'var') %no variable
                sepOnsetPeak(1:hd.nADCNumChannels, 1:hd.lActualEpisodes) = struct('r', []);
                %sepOP.r(1) - onset of SEP (ms from stimulus onset)
                %sepOP.r(2) - peak of SEP (ms from stimulus onset)
                %sepOP.r(3) - SEP amplitude (uV)
                %sepOP.r(4) - MUA during SEP
                %sepOP.r(5) - MUA after SEP
                %sepOP.r(6) - maximal rise rate (uV/ms, negative or positive)
                %sepOP.r(7) - average 3-point slope
            end
        end
        if treatAllCh %treat all channel
            fullChList = 1:(16 * floor(hd.nADCNumChannels / 16));%run over all LFP-channels
        elseif ~sopExist %no variable
            fullChList = rCh;%treat single requested channel
        else
            fullChList = [];%have not to treat SEPs
        end
        
        for rCh = fullChList %run over all channel (special calculation of maximal rise rate)
            for sw = 1:hd.lActualEpisodes %run over all sweeps
                sepOnsetPeak(rCh, sw).r = Inf(numel(zavp.realStim(sw).r), 7);%predefinition
            end
            segmEdge = [str2double(get(handles.SEPbegin, 'String')), str2double(get(handles.SEPend, 'String'))];%left and right shifts from synchro-point (ms)
            whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to synchro-events
            whtLFP = squeeze(whtLFP);%squeezed array
            
%             %common reference
%             if ~isempty((strfind(zavp.file, '2017-12-08')))
%                 jj = ismember(-100:250, segmEdge(1):segmEdge(2));
%                 whtLFP = whtLFP - repmat(wLFPcomm(jj, rCh), 1, size(whtLFP, 2));%common reference
%             end

            trcLen = size(whtLFP, 1);%time-length of signals (points)
        
            % === analyse average trace === %
            singlTrac = smooth(mean(whtLFP, 2), 3);%smoothing
            singlTrac = ZavMultiDiff(singlTrac, 2:4) ./ ([ones(1, 10), (1:(trcLen - 10))]' .^ 0.3);%multistep differential
            [~, midPnt] = max(abs(singlTrac));%point of maximal rise speed (approximate half-rise point)
            frontSign = sign(singlTrac(midPnt));%signature in the ponit of maximal slope
            jj = midPnt + (-5:5)'; jj = jj((jj > 0) & (jj <= trcLen));%index points of current trace to be treated
            singlTrac = mean(whtLFP, 2) * frontSign;%forced positive polarity
            ii = diff(singlTrac(jj));%differential (speed of voltage changes)
            jj01 = jj;%copy of index
            jj01(ii <= 0) = -1 + rand(sum(ii <= 0), 1);%riples for wrong sign
            jj01(ii > 0) = 1;%flate for right sign
            [pntsL, pntsR] = ZavFindFlats(jj01);%flat segments

            segmEdge(2) = Inf;
            if (~isempty(pntsL) && ~isempty(pntsR)) %at least one flat segment exist
                [~, z] = max(diff([pntsL, pntsR], [], 2));%longest segment of uninterrupted rise
                ii = ii(pntsL(z):pntsR(z));%only negative differences
                jj = (jj(pntsL(z)):(jj(pntsR(z)) + 1))';%indices of rise-front points
                midPnt = round(interp1(singlTrac(jj), jj, mean(singlTrac(jj([1, end])))));%median point %1)midPnt = round(mean(jj([pntsL(z), pntsR(z)])));%approximate half-rise point

                z = max(ii) / 20;%threshold for rise rate near onset
                p1 = find(diff(singlTrac(1:(midPnt + 0))) < z, 1, 'last');
                if ~isempty(p1)
                    p1 = segmEdge(1) + p1;%onset of SEP
                else
                    p1 = segmEdge(1);
                end

                z = max(ii) / 20;%threshold for rise rate near peak
                p2 = find(diff(singlTrac(midPnt:end)) < z, 1, 'first');
                if ~isempty(p2)
                    p2 = segmEdge(1) + p2 + midPnt - 2;%onset of SEP
                end
                
                %new segmEdge, riseThreshold, peakThreshold
                segmEdge = [p1 - 5, p2 + 5];%left and right shifts from synchro-point (ms)
                %segmEdge(1) = max(segmEdge(1), 1);
                %riseThreshold = [-max(ii) / 15, max(ii) / 1];%threshold for rise rate near onset
                %peakThreshold = [max(ii) / 8, -max(ii) / 1.5];%threshold for rise rate near onset
            end
            % === end of analysis of average trace === %
            
            if (all(isfinite(segmEdge)) && (numel(segmEdge) == 2))
                %new segments of LFP (within new range)
                whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to synchro-events
                whtLFP = squeeze(whtLFP);%squeezed array
                
%                 %common reference
%                 if ~isempty((strfind(zavp.file, '2017-12-08')))
%                     jj = ismember(-100:250, segmEdge(1):segmEdge(2));
%                     whtLFP = whtLFP - repmat(wLFPcomm(jj, rCh), 1, size(whtLFP, 2));%common reference
%                 end

                trcLen = size(whtLFP, 1);%time-length of signals (points)
                for segmIndx = 1:size(segms, 1) %run over segments %sepOnsetPeak(rCh, sw).r, 1)
                    sw = segms(segmIndx, 3);%current sweep (number of sweep where syncro-event was detected)
                    g = segms(segmIndx, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix

                    %singlTrac = smooth(whtLFP(:, segmIndx), 3);%smoothing
                    %singlTrac = ZavMultiDiff(singlTrac, 2:4) ./ ([ones(1, 10), (1:(trcLen - 10))]' .^ 0.3);%multistep differential
                    singlTrac = diff(whtLFP(:, segmIndx));%differential
                    [~, midPnt] = max(abs(singlTrac));%point of maximal rise speed (approximate half-rise point)
                    frontSign = sign(singlTrac(midPnt));%signature in the ponit of maximal slope
                    jj = midPnt + (-5:5)'; jj = jj((jj > 0) & (jj <= trcLen));%index points of current trace to be treated
                    ii = diff(frontSign * whtLFP(jj, segmIndx));%differential (speed of voltage changes, forced positive polarity)
                    jj01 = jj;%copy of index
                    jj01(ii <= 0) = -1.1 + rand(sum(ii <= 0), 1);%riples for wrong sign
                    jj01(ii > 0) = 1;%flate for right sign
                    [pntsL, pntsR] = ZavFindFlats(jj01);%flat segments

                    if (~isempty(pntsL) && ~isempty(pntsR)) %at least one flat segment exist
                        p2 = mean(ii(pntsL(1):pntsR(1))); p1 = 1;
                        for z = 2:length(pntsL) %run over flat segments
                            if (mean(ii(pntsL(z):pntsR(z))) > p2)
                                p2 = mean(ii(pntsL(z):pntsR(z)));
                                p1 = z;
                            end
                        end
                        z = p1;%steepest segment
                        %[~, z] = max(diff([pntsL, pntsR], [], 2));%longest segment of uninterrupted rise
                            
                        ii = ii(pntsL(z):pntsR(z));%sharp front only
                        jj = (jj(pntsL(z)):(jj(pntsR(z)) + 1))';%indices of rise-front points
                        %midPnt = round(interp1(whtLFP(jj, segmIndx), jj, mean(whtLFP(jj([1, end]), segmIndx))));%median point %1)midPnt = round(mean(jj([pntsL(z), pntsR(z)])));%approximate half-rise point
                        [~, midPnt] = max(ii);
                        midPnt = jj(midPnt);

                        singlTrac = frontSign * whtLFP((midPnt + 0):-1:1, segmIndx);%inverted time, forced positive polarity
                        riseThreshold = -max(ii) / 10;%threshold for rise rate (near onset)
                        p1 = find(diff(singlTrac) > riseThreshold, 1, 'first');%decrease of slope
                        %          ([0; diff(singlTrac, 2)] > riseThreshold(2)), ... %or big jump of rate (strong relative changes)
                        if ~isempty(p1)
                            sepOnsetPeak(rCh, sw).r(g, 1) = segmEdge(1) + (midPnt - p1);%onset of SEP
                        end
                        
                        singlTrac = frontSign * whtLFP(midPnt:end, segmIndx);%direct (right) time, forced positive polarity
                        peakThreshold = max(ii) / 3;%threshold for rise rate (near peak)
                        p2 = find(diff(singlTrac) < peakThreshold, 1, 'first');%absolute fall rate
                        %          ([0; diff(whtLFP(midPnt:end, segmIndx), 2)] < peakThreshold(2)), ... %or big jump of rate (strong relative changes)
                        if ~isempty(p2)
                            sepOnsetPeak(rCh, sw).r(g, 2) = segmEdge(1) + p2 + midPnt - 2;%peak of SEP
                        end

                        sepOnsetPeak(rCh, sw).r(g, 6) = max(ii) * frontSign;%maximal rise rate (uV/ms, negative or positive)
                        [~, z] = max(ii);
                        if (((jj(z) - 1) > 0) && ((jj(z) + 1) <= trcLen)) %
                            if ((jj(z) - 1) > (midPnt - p1))
                                sepOnsetPeak(rCh, sw).r(g, 7) = diff(whtLFP(jj(z) + [-1, 1], segmIndx)) / 2;%average 3-point slope
                            elseif ((jj(z) + 2) <= trcLen)
                                sepOnsetPeak(rCh, sw).r(g, 7) = diff(whtLFP(jj(z) + [0, 2], segmIndx)) / 2;%average 3-point slope
                            end
                        end
                    end

                    if all(isfinite(sepOnsetPeak(rCh, sw).r(g, 1:2))) % && all(sepOnsetPeak(rCh, sw).r(g, 1:2) > 0))
                        p1 = sepOnsetPeak(rCh, sw).r(g, 1) - segmEdge(1) + 1;%onset point
                        p2 = sepOnsetPeak(rCh, sw).r(g, 2) - segmEdge(1) + 1;%peak point
                        sepOnsetPeak(rCh, sw).r(g, 3) = diff(whtLFP([p1, p2], segmIndx));%amplitude
                        
                        %singlTrac = ZavSynchLFP(zavp, hd, segms(segmIndx, :), [0, 30], lfp, rCh, false);
                        %plot(0:30, singlTrac, '.-', sepOnsetPeak(rCh, sw).r(g, 1:2), singlTrac(sepOnsetPeak(rCh, sw).r(g, 1:2) + 1), 'or', jj(z) + segmEdge(1) - 1, singlTrac(jj(z) + segmEdge(1)), 'sg')
                        %disp([rCh, segmIndx])
                    end
                end
            end
        end
        if (treatAllCh || treatAllRecs) %all channels automatical treatment
            lfpVar = getappdata(expAn, 'lfpVar');
            spks = getappdata(expAn, 'spks');
            %sepOnsetPeak2 = sepOnsetPeak;%automatically calculated (maximal rise rate confirmed)
            sepOnsetPeak2 = GetSpikesCount(sepOnsetPeak, hd, spks, lfpVar, segms, handles);%calculate spikes density near SEP
            save([pth, flNm], 'sepOnsetPeak2', '-append')%save resultes
        else
            setappdata(expAn, 'sepOnsetPeak', sepOnsetPeak);%SEP parameters
        end
        
        if ~treatAllRecs %treat single record
            rCh = getappdata(expAn, 'rCh');%principal channel
            set(handles.AcceptSweep, 'Enable', 'on')
            set(handles.AcceptAll, 'Enable', 'on')
            set(handles.SkipSweep, 'Enable', 'on')
            set(handles.BackSweep, 'Enable', 'on')

            cla(ax1), hold(ax1, 'on')
    %         p1 = segmEdge(1) + 5; p2 = segmEdge(2) - 5;
            segmEdge = [str2double(get(handles.SEPbegin, 'String')), str2double(get(handles.SEPend, 'String'))];%left and right shifts from synchro-point (ms)
            setappdata(expAn, 'segmEdge', segmEdge);%left and right shifts from synchro-point (ms)
            segmEdge = segmEdge + [-5, 10];%left and right shifts from synchro-point (ms)
            whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to synchro-events
            whtLFP = squeeze(whtLFP);%squeezed array
            
%             %common reference
%             if ~isempty((strfind(zavp.file, '2017-12-08')))
%                 jj = ismember(-100:250, segmEdge(1):segmEdge(2));
%                 whtLFP = whtLFP - repmat(wLFPcomm(jj, rCh), 1, size(whtLFP, 2));%common reference
%             end

            setappdata(expAn, 'whtLFP', whtLFP);%LFP time course
            zavp.prm.Fs = 1 / (zavp.siS * zavp.rarStep);%sampling frequency for rare data
            tm = segmEdge(1):segmEdge(2);%time matrix for raw or downsampled data (ms)

    %         plot(ax1, tm, mean(whtLFP, 2), 'Color', 0.0 * ones(1, 3), 'LineWidth', 2)
    %         plot(ax1, tm, [mean(whtLFP, 2) - std(whtLFP, [], 2), mean(whtLFP, 2) + std(whtLFP, [], 2)], 'Color', 0.9 * ones(1, 3), 'LineWidth', 2)
    %         plot(ax1, repmat([p1, p2], 2, 1), repmat(ylim' * 2, 1, 2), 'm') 

            setappdata(expAn, 'tm', tm);%time vector
            segmIndx = 1;%segments counter
            set(textMsgs, 'String', ['segment ', num2str(segmIndx), ' of ', num2str(size(segms, 1))])%count of segments (string)
            sw = segms(segmIndx, 3);%number of sweep where syncro-event was detected
            g = segms(segmIndx, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
            p1 = sepOnsetPeak(rCh, sw).r(g, 1);%onset points
            p2 = sepOnsetPeak(rCh, sw).r(g, 2);%peak points
            [~, ii(1)] = min(abs(tm - p1)); [~, ii(2)] = min(abs(tm - p2));

            nulLevel = 0;%mean(whtLFP(1:(ii(1) - 2), segmIndx));%null-level
            tracH = plot(ax1, tm, whtLFP(:, segmIndx) - nulLevel, 'b.-');%plot first trace and create handle of the line
            set(tracH, 'ButtonDownFcn', {@Axes1_ButtonDownFcn, handles})%define reaction on mouse click
            setappdata(expAn, 'tracH', tracH);%handle of trace
            g = whtLFP(ii(1:2), segmIndx) - nulLevel;%principal points (onset and peak)
            if sopExist %variable exist
                pntsH = plot(ax1, [p1, p2], g, 'sg');%add onset-peak points and create handle of the graph
            else %first definition of variable
                pntsH = plot(ax1, [p1, p2], g, 'or');%add onset-peak points and create handle of the graph
            end
            setappdata(expAn, 'pntsH', pntsH);%handle of points
            setappdata(expAn, 'segmIndx', segmIndx);%segment index
            set(ax1, 'XLim', tm([1, end]));
            set(ax1, 'YLim', [min(whtLFP(:, segmIndx)) - 5, max(whtLFP(:, segmIndx)) + 5] - nulLevel)
%             if all(isfinite(g)) %onset and peak points are defined correctly
%                 set(ax1, 'YLim', sort(mean(g) + (diff(g) / 1) * [-1, 1]));%horizontal and vertical axes limits
%             end
        end
    else %no evoked signals
        msgs = get(textMsgs, 'String');
        msgs{end + 1} = 'no evoked signals';%display message
        set(textMsgs, 'String', msgs);
    end
    
    %OFF-response of evoked signals
    if ~isempty(zavp.realStim(1).f) %evoked signals (off-response)
        pth = getappdata(expAn, 'pth');%directory
        %hd = getappdata(expAn, 'hd');%header
        lfp = getappdata(expAn, 'lfp');%LFP
        zavp = getappdata(expAn, 'zavp');%metadata
        treatAllCh = isequal(get(handles.TreatAllChannels, 'State'), 'on');%treat all channel
        set(expAn, 'Name', ['Express analysis - ', pth, flNm]);%full name
        segms = GetSegments(hd, zavp, 2);%segments of evoked activity
        setappdata(expAn, 'segms', segms);%segments
        rCh = str2double(get(handles.PrincChan, 'String'));%zavp.chL4(1);%principal channel
        setappdata(expAn, 'rCh', rCh);%principal channel

        ii = [];%principal channels
        sepOnsetPeak3(1:hd.nADCNumChannels, 1:hd.lActualEpisodes) = struct('r', []);

        if treatAllCh %treat all channel
            fullChList = 1:(16 * floor(hd.nADCNumChannels / 16));%run over all LFP-channels
        elseif ~sopExist %no variable
            fullChList = rCh;%treat single requested channel
        else
            fullChList = [];%have not to treat SEPs
        end
        
        for rCh = fullChList %run over all channel (special calculation of maximal rise rate)
            for sw = 1:hd.lActualEpisodes %run over all sweeps
                sepOnsetPeak3(rCh, sw).r = Inf(numel(zavp.realStim(sw).f), 7);%predefinition
            end
            segmEdge = [str2double(get(handles.SEPbegin, 'String')), str2double(get(handles.SEPend, 'String'))];%left and right shifts from synchro-point (ms)
            whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to synchro-events
            whtLFP = squeeze(whtLFP);%squeezed array
            
%             %common reference
%             if ~isempty((strfind(zavp.file, '2017-12-08')))
%                 jj = ismember(-100:250, segmEdge(1):segmEdge(2));
%                 whtLFP = whtLFP - repmat(wLFPcomm(jj, rCh), 1, size(whtLFP, 2));%common reference
%             end

            trcLen = size(whtLFP, 1);%time-length of signals (points)
        
            % === analyse average trace === %
            singlTrac = smooth(mean(whtLFP, 2), 3);%smoothing
            singlTrac = ZavMultiDiff(singlTrac, 2:4) ./ ([ones(1, 10), (1:(trcLen - 10))]' .^ 0.3);%multistep differential
            [~, midPnt] = max(abs(singlTrac));%point of maximal rise speed (approximate half-rise point)
            frontSign = sign(singlTrac(midPnt));%signature in the ponit of maximal slope
            jj = midPnt + (-5:5)'; jj = jj((jj > 0) & (jj <= trcLen));%index points of current trace to be treated
            singlTrac = mean(whtLFP, 2) * frontSign;%forced positive polarity
            ii = diff(singlTrac(jj));%differential (speed of voltage changes)
            jj01 = jj;%copy of index
            jj01(ii <= 0) = -1 + rand(sum(ii <= 0), 1);%riples for wrong sign
            jj01(ii > 0) = 1;%flate for right sign
            [pntsL, pntsR] = ZavFindFlats(jj01);%flat segments

            segmEdge(2) = Inf;
            if (~isempty(pntsL) && ~isempty(pntsR)) %at least one flat segment exist
                [~, z] = max(diff([pntsL, pntsR], [], 2));%longest segment of uninterrupted rise
                ii = ii(pntsL(z):pntsR(z));%only negative differences
                jj = (jj(pntsL(z)):(jj(pntsR(z)) + 1))';%indices of rise-front points
                midPnt = round(interp1(singlTrac(jj), jj, mean(singlTrac(jj([1, end])))));%median point %1)midPnt = round(mean(jj([pntsL(z), pntsR(z)])));%approximate half-rise point

                z = max(ii) / 20;%threshold for rise rate near onset
                p1 = find(diff(singlTrac(1:(midPnt + 0))) < z, 1, 'last');
                if ~isempty(p1)
                    p1 = segmEdge(1) + p1;%onset of SEP
                else
                    p1 = segmEdge(1);
                end

                z = max(ii) / 20;%threshold for rise rate near peak
                p2 = find(diff(singlTrac(midPnt:end)) < z, 1, 'first');
                if ~isempty(p2)
                    p2 = segmEdge(1) + p2 + midPnt - 2;%onset of SEP
                end
                
                %new segmEdge, riseThreshold, peakThreshold
                segmEdge = [p1 - 5, p2 + 5];%left and right shifts from synchro-point (ms)
                segmEdge(1) = max(segmEdge(1), 1);
                %riseThreshold = [min(ii) / 15, -min(ii) / 1];%threshold for rise rate near onset
                %peakThreshold = [min(ii) / 8, -min(ii) / 1.5];%threshold for rise rate near onset
            end
            % === end of analysis of average trace === %
            
            if (all(isfinite(segmEdge)) && (numel(segmEdge) == 2))
                %new segments of LFP (within new range)
                whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to synchro-events
                whtLFP = squeeze(whtLFP);%squeezed array
                
%                 %common reference
%                 if ~isempty((strfind(zavp.file, '2017-12-08')))
%                     jj = ismember(-100:250, segmEdge(1):segmEdge(2));
%                     whtLFP = whtLFP - repmat(wLFPcomm(jj, rCh), 1, size(whtLFP, 2));%common reference
%                 end

                trcLen = size(whtLFP, 1);%time-length of signals (points)
                for segmIndx = 1:size(segms, 1) %run over segments %sepOnsetPeak3(rCh, sw).r, 1)
                    sw = segms(segmIndx, 3);%current sweep (number of sweep where syncro-event was detected)
                    g = segms(segmIndx, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix

                    %singlTrac = smooth(whtLFP(:, segmIndx), 3);%smoothing
                    %singlTrac = ZavMultiDiff(singlTrac, 2:4) ./ ([ones(1, 10), (1:(trcLen - 10))]' .^ 0.3);%multistep differential
                    singlTrac = diff(whtLFP(:, segmIndx));%differential
                    [~, midPnt] = max(abs(singlTrac));%point of maximal rise speed (approximate half-rise point)
                    frontSign = sign(singlTrac(midPnt));%signature in the ponit of maximal slope
                    jj = midPnt + (-5:5)'; jj = jj((jj > 0) & (jj <= trcLen));%index points of current trace to be treated
                    ii = diff(frontSign * whtLFP(jj, segmIndx));%differential (speed of voltage changes, forced positive polarity)
                    jj01 = jj;%copy of index
                    jj01(ii <= 0) = -1 + rand(sum(ii <= 0), 1);%riples for wrong sign
                    jj01(ii > 0) = 1;%flate for right sign
                    [pntsL, pntsR] = ZavFindFlats(jj01);%flat segments

                    if (~isempty(pntsL) && ~isempty(pntsR)) %at least one flat segment exist
                        [~, z] = max(diff([pntsL, pntsR], [], 2));%longest segment of uninterrupted rise
                        ii = ii(pntsL(z):pntsR(z));%only negative differences
                        jj = (jj(pntsL(z)):(jj(pntsR(z)) + 1))';%indices of rise-front points
                        %midPnt = round(interp1(whtLFP(jj, segmIndx), jj, mean(whtLFP(jj([1, end]), segmIndx))));%median point %1)midPnt = round(mean(jj([pntsL(z), pntsR(z)])));%approximate half-rise point
                        [~, midPnt] = max(ii);
                        midPnt = jj(midPnt);

                        singlTrac = frontSign * whtLFP((midPnt + 0):-1:1, segmIndx);%inverted time, forced positive polarity
                        riseThreshold = -max(ii) / 10;%threshold for rise rate (near onset)
                        p1 = find(diff(singlTrac) > riseThreshold, 1, 'first');%decrease of slope
                        %          ([0; diff(singlTrac, 2)] > riseThreshold(2)), ... %or big jump of rate (strong relative changes)
                        if ~isempty(p1)
                            sepOnsetPeak3(rCh, sw).r(g, 1) = segmEdge(1) + (midPnt - p1);%onset of SEP
                        end

                        singlTrac = frontSign * whtLFP(midPnt:end, segmIndx);%direct (right) time, forced positive polarity
                        peakThreshold = max(ii) / 3;%threshold for rise rate (near peak)
                        p2 = find(diff(singlTrac) < peakThreshold, 1, 'first');%absolute rate fall
                        %          ([0; diff(whtLFP(midPnt:end, segmIndx), 2)] < peakThreshold(2)), ... %or big jump of rate (strong relative changes)
                        if ~isempty(p2)
                            sepOnsetPeak3(rCh, sw).r(g, 2) = segmEdge(1) + p2 + midPnt - 2;%peak of SEP
                        end

                        sepOnsetPeak3(rCh, sw).r(g, 6) = max(ii) * frontSign;%maximal rise rate (uV/ms, negative or positive)
                        [~, z] = max(ii);
                        if (((jj(z) - 1) > 0) && ((jj(z) + 1) <= trcLen)) %
                            if ((jj(z) - 1) > (midPnt - p1))
                                sepOnsetPeak3(rCh, sw).r(g, 7) = diff(whtLFP(jj(z) + [-1, 1], segmIndx)) / 2;%average 3-point slope
                            elseif ((jj(z) + 2) <= trcLen)
                                sepOnsetPeak3(rCh, sw).r(g, 7) = diff(whtLFP(jj(z) + [0, 2], segmIndx)) / 2;%average 3-point slope
                            end
                        end
                    end

                    if all(isfinite(sepOnsetPeak3(rCh, sw).r(g, 1:2))) % && all(sepOnsetPeak3(rCh, sw).r(g, 1:2) > 0))
                        p1 = sepOnsetPeak3(rCh, sw).r(g, 1) - segmEdge(1) + 1;%onset point
                        p2 = sepOnsetPeak3(rCh, sw).r(g, 2) - segmEdge(1) + 1;%peak point
                        sepOnsetPeak3(rCh, sw).r(g, 3) = diff(whtLFP([p1, p2], segmIndx));%amplitude
                        
                        %singlTrac = ZavSynchLFP(zavp, hd, segms(segmIndx, :), [0, 30], lfp, rCh, false);
                        %plot(0:30, singlTrac, '.-', sepOnsetPeak3(rCh, sw).r(g, 1:2), singlTrac(sepOnsetPeak3(rCh, sw).r(g, 1:2) + 1), 'or', jj(z) + segmEdge(1) - 1, singlTrac(jj(z) + segmEdge(1)), 'sg')
                        %disp([rCh, segmIndx])
                    end
                end
            end
        end
        lfpVar = getappdata(expAn, 'lfpVar');
        spks = getappdata(expAn, 'spks');
        sepOnsetPeak3 = GetSpikesCount(sepOnsetPeak3, hd, spks, lfpVar, segms, handles);%calculate spikes density near SEP
        save([pth, flNm], 'sepOnsetPeak3', '-append')%save resultes
    else %no evoked signals
        msgs = get(textMsgs, 'String');
        msgs{end + 1} = 'no evoked signals';%display message
        set(textMsgs, 'String', msgs);
    end
end



function [sepOnsetPeak, ii] = CorrectSepOP(expAn, sepOnsetPeak, hd, zavp)
%correct sepOnsetPeak
ii = [];%principal channels
if (numel(sepOnsetPeak) == (hd.nADCNumChannels * hd.lActualEpisodes)) %all-channels-all-sweeps structure
    ii = [];%processed principal channels
    for t = 1:hd.nADCNumChannels %run over channels
        if ~isempty(vertcat(sepOnsetPeak(t, :).r))
            ii = [ii, t];%processed principal channel
        end
    end
elseif (numel(sepOnsetPeak) == 1) %single-channel-single-sweep structure
    ii = sepOnsetPeak(1, 1).r;%copy of SEP parameters
    sepOnsetPeak(1:hd.nADCNumChannels, 1:hd.lActualEpisodes) = struct('r', []);%convert to all-channels-single-sweep structure
    sepOnsetPeak(zavp.chL4(1), 1).r = ii;%convert to all-channels-single-sweep structure
    ii = zavp.chL4(1);%principal channel
    setappdata(expAn, 'sepOnsetPeak', sepOnsetPeak);%earlier calculated SEP parameters (reshaped)
else %notall-channels-notall-sweeps structure
    errMsgH = errordlg('Ask Andrey Zakharov to fix the program', 'Unexpected data format', 'modal');
    return %no conversion provided
end



function sepOnsetPeak = GetSpikesCount(sepOnsetPeak, hd, spks, lfpVar, segms, handles)
%calculate spikes density near SEP
fullChList = 1:(16 * floor(hd.nADCNumChannels / 16));%run over all LFP-channels %rCh = getappdata(expAn, 'rCh');%principal channel
for rCh = fullChList %run over all LFP-channels
    if ~isempty(vertcat(sepOnsetPeak(rCh, :).r)) %SEP parameters calculated
        prgMult = (-1) * str2double(get(handles.Porog, 'String'));%threshold multiplier
        for sw = 1:hd.lActualEpisodes %run over sweeps
           jj = spks(rCh, sw).ampl <= (prgMult * min(lfpVar(rCh, :)));
           spks(rCh, sw).tStamp = spks(rCh, sw).tStamp(jj);
           %spks(rCh, sw).ampl = spks(rCh, sw).ampl(jj);
           %spks(rCh, sw).shape = spks(rCh, sw).shape(:, jj);
        end
        sgmNum = size(segms, 1);%number of segments
        spkDens = Inf(sgmNum, 2);%MUA density
        for t = 1:sgmNum %run over segments (stimulus)
            sw = segms(t, 3);%number of sweep
            g = segms(t, 4);%number of stimulus (in matrix zavp.realStim)
            if (~isempty(spks(rCh, sw).tStamp) && all(isfinite(sepOnsetPeak(rCh, sw).r(g, 1:3))))
                segmEdge = segms(t, 1) + sepOnsetPeak(rCh, sw).r(g, 1) + [0, 20];%SEP period
                spkDens(t, 1) = sum((spks(rCh, sw).tStamp >= segmEdge(1)) & (spks(rCh, sw).tStamp <= segmEdge(2)));%mua during SEP
                segmEdge = segms(t, 1) + sepOnsetPeak(rCh, sw).r(g, 1) + [20, 520];%period after SEP (burst)
                spkDens(t, 2) = sum((spks(rCh, sw).tStamp >= segmEdge(1)) & (spks(rCh, sw).tStamp <= segmEdge(2)));%mua after SEP (burst)
                sepOnsetPeak(rCh, sw).r(g, 4) = spkDens(t, 1);%MUA during SEP
                sepOnsetPeak(rCh, sw).r(g, 5) = spkDens(t, 2);%MUA after SEP
            end
        end
    end
end


% --- Executes on button press in AcceptSweep.
function AcceptSkipSweep_Callback(hObject, eventdata, handles)
% hObject    handle to AcceptSweep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%accept SEP parameters on current sweep
expAn = handles.ExpressAnalysis;%main form handle
ax1 = handles.Axes1;%handle of axes
textMsgs = handles.TextMsgs;%messeges handle
zavp = getappdata(expAn, 'zavp');%file accessories
whtLFP = getappdata(expAn, 'whtLFP');%LFP segmented
segms = getappdata(expAn, 'segms');%segments
segmIndx = getappdata(expAn, 'segmIndx');%segment index
sepOnsetPeak = getappdata(expAn, 'sepOnsetPeak');%SEP parameters
rCh = getappdata(expAn, 'rCh');%principal channel
tracH = getappdata(expAn, 'tracH');%handle of trace
pntsH = getappdata(expAn, 'pntsH');%handle of onset-peak points
tm = getappdata(expAn, 'tm');%time vector

sw = segms(segmIndx, 3);%number of sweep
g = segms(segmIndx, 4);%stimulus number
if isequal(get(hObject, 'Tag'), 'AcceptSweep') %accept sweep
    p = sort(get(pntsH, 'XData'));%key points coordinates
    sepOnsetPeak(rCh, sw).r(g, 1:2) = p;
    [~, ii(1)] = min(abs(tm - p(1)));
    [~, ii(2)] = min(abs(tm - p(2)));
    sepOnsetPeak(rCh, sw).r(g, 3) = diff(whtLFP(ii(1:2), segmIndx));
    segmIndx = segmIndx + 1;%next segment
elseif isequal(get(hObject, 'Tag'), 'SkipSweep') %reject sweep
    sepOnsetPeak(rCh, sw).r(g, :) = Inf(1, 6);%incalculable signal
    segmIndx = segmIndx + 1;%next segment
else %if isequal(get(hObject, 'Tag'), 'BackSweep') %one sweep back
    segmIndx = segmIndx - 1;%previous segment
    if (segmIndx < 1)
        segmIndx = 1;
    end
end
setappdata(expAn, 'sepOnsetPeak', sepOnsetPeak);%manually corrected SEP parameters

%draw next sweep
if (segmIndx > length(segms)) %all sweeps was treated
    lfpVar = getappdata(expAn, 'lfpVar');%
    spks = getappdata(expAn, 'spks');
    hd = getappdata(expAn, 'hd');
    sepOnsetPeak = GetSpikesCount(sepOnsetPeak, hd, spks, lfpVar, segms, handles);%recalculate spikes density near SEP
    pth = getappdata(expAn, 'pth');%directory
    flNm = getappdata(expAn, 'flNm');%filename
    save([pth, flNm], 'sepOnsetPeak', '-append')%save resultes
    
    cla(ax1)
    set(textMsgs, 'String', {});
    %set(expAn, 'Name', 'Express analysis');%erase full name
    return;
end
sw = segms(segmIndx, 3);%number of sweep
g = segms(segmIndx, 4);%stimulus number
p = sepOnsetPeak(rCh, sw).r(g, 1:2);%next onset and peak points
[~, ii(1)] = min(abs(tm - p(1)));%corresponding point of time vector
[~, ii(2)] = min(abs(tm - p(2)));%corresponding point of time vector
nulLevel = 0;%mean(whtLFP(1:(ii(1) - 2), segmIndx));%null-level
set(tracH, 'XData', tm, 'YData', whtLFP(:, segmIndx) - nulLevel);
g = whtLFP(ii(1:2), segmIndx) - nulLevel;%principal points (onset and peak)
set(pntsH, 'XData', tm(ii(1:2)), 'YData', g);%onset-peak points
set(ax1, 'XLim', tm([1, end]));
set(ax1, 'YLim', [min(whtLFP(:, segmIndx)) - 5, max(whtLFP(:, segmIndx)) + 5] - nulLevel)
% if (all(isfinite(g)) && (g(1) ~= g(2))) %onset and peak points are defined correctly
%     set(ax1, 'YLim', sort(mean(g) + max((abs(diff(g)) / 1), 150) * [-1, 1]))
% end
set(textMsgs, 'String', ['segment ', num2str(segmIndx), ' of ', num2str(size(segms, 1))]);
setappdata(expAn, 'segmIndx', segmIndx);%segment index



% --- Executes on button press in AcceptAll.
function AcceptAll_Callback(hObject, eventdata, handles)
% hObject    handle to AcceptAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

expAn = handles.ExpressAnalysis;%main form handle
segms = getappdata(expAn, 'segms');%segments
sepOnsetPeak = getappdata(expAn, 'sepOnsetPeak');%SEP parameters
whtLFP = getappdata(expAn, 'whtLFP');%LFP segmented
tm = getappdata(expAn, 'tm');%time vector
rCh = getappdata(expAn, 'rCh');%principal channel
tracH = getappdata(expAn, 'tracH');%handle of trace
pntsH = getappdata(expAn, 'pntsH');%handle of onset-peak points

segmIndx = length(segms);%go final sweep and accept
setappdata(expAn, 'segmIndx', segmIndx);%segment index
sw = segms(segmIndx, 3);%number of sweep
g = segms(segmIndx, 4);%stimulus number
p = sepOnsetPeak(rCh, sw).r(g, 1:2);%next onset and peak points
[~, ii(1)] = min(abs(tm - p(1)));%corresponding point of time vector
[~, ii(2)] = min(abs(tm - p(2)));%corresponding point of time vector
nulLevel = 0;%mean(whtLFP(1:(ii(1) - 2), segmIndx));%null-level
set(tracH, 'XData', tm, 'YData', whtLFP(:, segmIndx) - nulLevel);
g = whtLFP(ii(1:2), segmIndx) - nulLevel;%principal points (onset and peak)
set(pntsH, 'XData', tm(ii(1:2)), 'YData', g);%onset-peak points
AcceptSkipSweep_Callback(handles.AcceptSweep, [], handles)%call "Accept" button


% --- Executes on mouse press over axes background.
function Axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to Axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%reaction of mouse button down
expAn = handles.ExpressAnalysis;%main form handle
whtLFP = getappdata(expAn, 'whtLFP');%LFP segmented
segmIndx = getappdata(expAn, 'segmIndx');%segment index
ax1 = handles.Axes1;%handle of axes
tm = getappdata(expAn, 'tm');%time vector
pntsH = getappdata(expAn, 'pntsH');
p = get(pntsH, 'XData');%key points coordinates

xy = get(ax1, 'CurrentPoint');
x_lim = get(ax1, 'XLim'); y_lim = get(ax1, 'YLim');
clkInAx = ((xy(1, 1) >= x_lim(1)) && (xy(1, 1) <= x_lim(2)) && ...
           (xy(1, 2) >= y_lim(1)) && (xy(1, 2) <= y_lim(2)));
if (clkInAx)
    if isequal(get(expAn, 'SelectionType'), 'normal')
        z = 1;
    else
        z = 2;
    end
    p(z) = xy(1, 1);
    [~, ii(1)] = min(abs(tm - p(1))); [~, ii(2)] = min(abs(tm - p(2)));
    nulLevel = 0;%mean(whtLFP(1:(ii(1) - 2), segmIndx));%null-level
    set(pntsH, 'XData', p, 'YData', whtLFP(ii(1:2), segmIndx) - nulLevel)
end



function PrincChan_Callback(hObject, eventdata, handles)
% hObject    handle to PrincChan (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

expAn = handles.ExpressAnalysis;%main form handle
prncChan = handles.PrincChan;%principal channel handle
rCh = getappdata(expAn, 'rCh');%previous principal channel
new_rCh = str2num(get(prncChan, 'String'));%new principal channel
if isempty(new_rCh)
    set(prncChan, 'String', num2str(rCh));
    return;
end
hd = getappdata(expAn, 'hd');%header
if (new_rCh < 1) %wrong channel
    new_rCh = 1;%lowest channel
elseif ~isempty(hd) %header exist    
    if (new_rCh > hd.nADCNumChannels) %wrong channel
        new_rCh = num2str(hd.nADCNumChannels);%highest channel
    end
end
set(prncChan, 'String', new_rCh);
setappdata(expAn, 'rCh', new_rCh);%new principal channel
if ~isequal(new_rCh, rCh) %new channel
    AnalyseSinglEvoked(handles);%start analysis again
end



function SEP_BeginEnd_Callback(hObject, eventdata, handles)
% hObject    handle to SEPbegin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%change segment edges
expAn = handles.ExpressAnalysis;%main form handle
segmEdge = getappdata(expAn, 'segmEdge');%left and right shifts from synchro-point (ms)
if ~isequal(segmEdge, [str2double(get(handles.SEPbegin, 'String')), str2double(get(handles.SEPend, 'String'))]) %new segment edges
    rCh = getappdata(expAn, 'rCh');%previous principal channel
    sepOnsetPeak = getappdata(expAn, 'sepOnsetPeak');%SEP parameters
    for sw = 1:size(sepOnsetPeak, 2)
        sepOnsetPeak(rCh, sw).r = [];%
    end
    setappdata(expAn, 'sepOnsetPeak', sepOnsetPeak);%new SEP parameters
    AnalyseSinglEvoked(handles);%start analysis again
end



% --------------------------------------------------------------------
function ExportData_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%creat dialog for entering application time
if isempty(get(handles.PrincChan, 'String')) %request principal channel
    errMsgH = errordlg('Enter single channel as principal and try again', 'No channel specified', 'modal');
    return;%do not analysis
end
expAn = handles.ExpressAnalysis;%main form handle
nulTm = get(handles.ExportData, 'UserData');%null time (injection)
if isempty(nulTm) %no null time
    str1 = ''; str2 = ''; str3 = '';
else
    str1 = num2str(nulTm(1));%null-time (hours)
    str2 = num2str(nulTm(2));%null-time (minutes)
    str3 = num2str(nulTm(3));%number of "etalon" record
end
nulTimeFig = figure('Name', 'Enter injection time', 'numbertitle', 'off', 'Menubar', 'none', ...
    'Units', 'pixels', 'Position', [200 500 262 134], 'resize', 'off', 'ToolBar', 'none', ...
    'DoubleBuffer', 'on', 'Visible', 'off');
clr = get(nulTimeFig, 'Color');
hourEdtH = uicontrol('Parent', nulTimeFig, 'Style', 'edit', 'String', str1, 'Tag', 'Hour', 'FontSize', 12, ...
    'Units', 'normalized', 'Position', [0.339694656488549 0.619402985074627 0.263358778625954 0.246268656716418]);%, ...
minEdtH = uicontrol('Parent', nulTimeFig, 'Style', 'edit', 'String', str2, 'Tag', 'Minute', 'FontSize', 12, ...
    'Units', 'normalized', 'Position', [0.66030534351145 0.619402985074627 0.263358778625954 0.246268656716418]);%, ...
txtH1 = uicontrol('Parent', nulTimeFig, 'Style','text', 'Units', 'normalized', 'FontSize', 9,...
    'HorizontalAlignment', 'right', 'Position', [0.0534351145038168 0.634328358208955 0.25 0.223880597014925], ...
    'String', {'application'; 'time'}, 'Tag', 'text1', 'BackgroundColor', clr);
txtH2 = uicontrol('Parent', nulTimeFig, 'Style', 'text', 'String', ':', 'FontSize', 14, 'Tag','text2', ...
    'Units', 'normalized', 'Position', [0.610687 0.649253731343284 0.0299999999999999 0.186567164179104], ...
    'BackgroundColor', clr);
txtH3 = uicontrol('Parent', nulTimeFig, 'Style', 'text', 'Units', 'normalized', 'FontSize', 9,...
    'HorizontalAlignment', 'right', 'Position',[0.393129770992366 0.865671641791045 0.16412213740458 0.119402985074627], ...
    'String', 'hours', 'Tag', 'text3', 'BackgroundColor', clr);
txtH4 = uicontrol('Parent', nulTimeFig, 'Style','text', 'Units', 'normalized', 'FontSize', 9,...
    'HorizontalAlignment', 'right', 'Position',[0.66030534351145 0.865671641791045 0.194656488549618 0.119402985074627], ...
    'String', 'minutes', 'Tag','text4', 'BackgroundColor', clr);
etalonEdtH = uicontrol('Parent', nulTimeFig, 'Style', 'edit', 'String', str3, 'Tag', 'EtalonRec', 'FontSize', 12, ...
    'Units', 'normalized', 'Position', [0.49618320610687 0.350746268656716 0.209923664122137 0.194029850746269]);%, ...
txtH5 = uicontrol('Parent', nulTimeFig, 'Style', 'text', 'Units', 'normalized', 'FontSize', 9,...
    'HorizontalAlignment', 'right', 'Position', [0.209923664122137 0.335820895522388 0.251908396946565 0.223880597014925], ...
    'String', {'etalon'; 'record #'}, 'Tag', 'text5', 'BackgroundColor', clr);
plotGraphs = uicontrol('Parent', nulTimeFig, 'Style', 'pushbutton', 'String', 'Plot', 'Tag', 'PlotGraphs', 'FontSize', 9, ...
    'Units', 'normalized', 'Position', [0.687022900763358 0.0597014925373134 0.270992366412214 0.201492537313433], ...
    'Callback', {@CreatDynamicFigs, handles});
exportData = uicontrol('Parent', nulTimeFig, 'Style', 'pushbutton', 'String', 'Export (xls)', 'Tag', 'ExportData', 'FontSize', 9, ...
    'Units', 'normalized', 'Position', [0.049618320610687 0.0597014925373134 0.270992366412214 0.201492537313433], ...
    'Callback', {@CreatXLSfile, handles});

setappdata(expAn, 'hourEdtH', hourEdtH);%handle of hour-edit
setappdata(expAn, 'minEdtH', minEdtH);%handle of minute-edit
setappdata(expAn, 'etalonEdtH', etalonEdtH);%handle of etalon-record-edit
setappdata(expAn, 'nulTimeFig', nulTimeFig);%handle of dialog window
set(nulTimeFig , 'Visible', 'on')%show dialog

% function WriteNullTime(hObject, eventdata, handles)
% %send null time (application time) to main form
% hm = get(handles.ExportData, 'UserData');
% if isequal(get(hObject, 'Tag'), 'Hour') %hour
%     hm(1) = num2str(hObject, 'String');%entered hour
% elseif isequal(get(hObject, 'Tag'), 'Minute') %hour %minute
%     hm(2) = num2str(hObject, 'String');%entered minute
% end
% set(handles.ExportData, 'UserData', hm);%set time



function segms = GetSegments(hd, zavp, onoff)
%segments of evoked activity
segms = zeros(hd.lActualEpisodes * numel(vertcat(zavp.realStim(:).r)), 4);%preallocation of memory for stimuli moments
z = 1;%through counter of stimuli
if (onoff == 1) %on
    for sw = 1:hd.lActualEpisodes %run over sweeps
        for t = 1:numel(zavp.realStim(sw).r) %run over wanted synchro-events
            segms(z, 1) = zavp.realStim(sw).r(t) / zavp.rarStep;%position of synchro-event (ms from record begin (from 0))
            %segms(z, 1) = segms(z, 1);% + sepOnsetPeak(rCh, sw).r(z, 1);%shift to SEP onset
            segms(z, 2) = zavp.stimCh;%number of channel where syncro-event was detected
            segms(z, 3) = sw;%number of sweep where syncro-event was detected
            segms(z, 4) = t;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
            z = z + 1;%through counter of stimuli
        end
    end
else %off
    for sw = 1:hd.lActualEpisodes %run over sweeps
        for t = 1:numel(zavp.realStim(sw).f) %run over wanted synchro-events
            segms(z, 1) = zavp.realStim(sw).f(t) / zavp.rarStep;%position of synchro-event (ms from record begin (from 0))
            %segms(z, 1) = segms(z, 1);% + sepOnsetPeak(rCh, sw).r(z, 1);%shift to SEP onset
            segms(z, 2) = zavp.stimCh;%number of channel where syncro-event was detected
            segms(z, 3) = sw;%number of sweep where syncro-event was detected
            segms(z, 4) = t;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
            z = z + 1;%through counter of stimuli
        end
    end
end
segms(z:end, :) = [];%delete excess



% --------------------------------------------------------------------
function CreatXLSfile(hObject, eventdata, handles)
% hObject    handle to ExportData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%export data to xls-file
expAn = handles.ExpressAnalysis;%main form handle

hourEdtH = getappdata(expAn, 'hourEdtH');%handle of hour-edit
minEdtH = getappdata(expAn, 'minEdtH');%handle of minute-edit
etalonEdtH = getappdata(expAn, 'etalonEdtH');%handle of etalon-record-edit
nulTimeFig = getappdata(expAn, 'nulTimeFig');%handle of dialog window
nulTm = [str2double(get(hourEdtH, 'String')), str2double(get(minEdtH, 'String')), str2double(get(etalonEdtH, 'String'))];
rmappdata(expAn, 'hourEdtH');%remove handle
rmappdata(expAn, 'minEdtH');%remover handle
rmappdata(expAn, 'etalonEdtH');%remover handle
close(nulTimeFig)%close dialog

if ((numel(nulTm) < 3) || any(isnan(nulTm))) %no null time asigned
    return;
end
t0 = nulTm(1) * 3600 + nulTm(2) * 60;%injection time (minutes from day begin)

set(handles.ExportData, 'UserData', nulTm);%null time (injection)
fBase = getappdata(expAn, 'fBase');%get list of recordations
rList = get(handles.recList, 'Data');%list of wanted files
rList = vertcat(rList{:, 2});%flags of accepted files
load(fBase(1).matFile, 'hd')
fullChList = 1:(16 * floor(hd.nADCNumChannels / 16));%run over all LFP-channels
prgMult = (-1) * str2double(get(handles.Porog, 'String'));%threshold multiplier

load(fBase(nulTm(3)).matFile, 'lfpVar')
etPrgs = lfpVar;%etalon thresholds

k = find(fBase(1).matFile == '\', 1, 'last');
z = find(fBase(1).matFile(1:(k - 2)) == '\', 1, 'last');
pf = fBase(1).matFile(1:k);%pathname for xls saving
flNm = fBase(1).matFile((z + 1):(k - 1));%directory name

for rCh = fullChList %run over all channels
    paramToXLS = cell(5e3, 12);%memory preallocation
    paramToXLS(1, :) = {'dir', 'file', 'number', 'time of day', 'time after appl.', ...
              'SEP onset', 'SEP peak', 'SEP ampl', '0-20ms MUA', '20-520ms MUA', 'SEP slope', 'spnt MUA'};
    paramToXLS(2, :) = {'', '', '#', 'hh:mm:ss', 's', ...
              'ms', 'ms', 'uV', 'units', 'units', 'uV/ms', 'unit/s'};
    px = 3;
    for w = 1:length(fBase) %run over records
        if rList(w) %treat wanted files only
            k = find(fBase(w).matFile(1:(end - 2)) == '\', 1, 'last');
            varInFl = who(matfile(fBase(w).matFile));
            if any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
                sw = 1;
                if ismember('zp', varInFl)
                    load(fBase(w).matFile, 'zp')
                    zavp = zp;
                else
                    load(fBase(w).matFile, 'zavp')
                end
                if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
                    load(fBase(w).matFile, 'hd', 'sepOnsetPeak', 'sepOnsetPeak2', 'spks')
                    if ~isempty(sepOnsetPeak(rCh, sw).r)
                        sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
                    end
                elseif ismember('sepOnsetPeak', varInFl)
                    load(fBase(w).matFile, 'hd', 'sepOnsetPeak', 'spks')
                    sepOnsetPeak2 = sepOnsetPeak;
                elseif ismember('sepOnsetPeak2', varInFl)
                    load(fBase(w).matFile, 'hd', 'sepOnsetPeak2', 'spks')
                end
                segms = GetSegments(hd, zavp, 1);%segments (evoked)

                timeStr = repmat(' ', size(segms, 1), 8);%"Moscow" time of the segment begin (part of day)
                for z = 1:size(segms, 1) %run over segments
                    timeStr(z, :) = MoscowTime(hd.recTime(1) + (segms(z, 1) * 1e-3));%"Moscow" time of the segment begin (part of day)
                end

                sepOnsetPeak2 = GetSpikesCount(sepOnsetPeak2, hd, spks, etPrgs, segms, handles);%new spikes count with etalon threshold and wanted level
                segmN = size(segms, 1);%numbler of segments (stimuli)
                for z = 1:segmN %run over stimuli
                    sepOnsetPeak2(rCh, sw).r(z, ~isfinite(sepOnsetPeak2(rCh, sw).r(z, :))) = NaN;
                    if (z > 1)
                        spntMUA = sum((spks(rCh, sw).tStamp >= (segms(z - 1, 1) + 1000)) & (spks(rCh, sw).tStamp <= (segms(z, 1) - 0)));
                        spntMUA = spntMUA * 1e3 / ((segms(z, 1) - 0) - (segms(z - 1, 1) + 1000));%unit/s
                    else %(z == 1) %first stimulus
                        spntMUA = sum((spks(rCh, sw).tStamp >= 0) & (spks(rCh, sw).tStamp <= (segms(z, 1) - 0)));
                        spntMUA = spntMUA * 1e3 / (segms(z, 1) - 0);%unit/s
                    end
                    paramToXLS(px, :) = {fBase(w).matFile(1:k), fBase(w).matFile((k + 1):end), z, timeStr(z, :), ... (dir, file, number)
                        round((hd.recTime(1) - t0 + ((zavp.realStim(sw).r(z) - 1) / zavp.rarStep) * 1e-3) * 100) / 100, ... (time after appl.)
                        sepOnsetPeak2(rCh, sw).r(z, 1), sepOnsetPeak2(rCh, sw).r(z, 2), round(sepOnsetPeak2(rCh, sw).r(z, 3) * 10) / 10, ... (SEP onset, SEP peak, SEP ampl)
                        sepOnsetPeak2(rCh, sw).r(z, 4), sepOnsetPeak2(rCh, sw).r(z, 5), round(sepOnsetPeak2(rCh, sw).r(z, 6) * 10) / 10, ... (0-20ms MUA, 20-520ms MUA, SEP slope)
                        round(spntMUA * 100) / 100 ... spontaneous MUA
                        };
                    px = px + 1;
                end
            end
        end
    end
    paramToXLS(px:end, :) = [];%delete excess
    xlswrite([pf, flNm, '.xlsx'], paramToXLS, ['ch_', num2str(rCh)]);%export data to xls-file
end



function timeStr = MoscowTime(timeNum)
%convert time as number to string in format hh:mm
%INPUTS
%timeNum - time as number (seconds from day beginning)
%OUTPUTS
%timeStr - time as string in format hh:mm

tm = zeros(1, 2);
tm(1) = floor(timeNum / 3600);%hours
tm(2) = floor((timeNum - (tm(1) * 3600)) / 60);%minutes
tm(3) = floor(timeNum - (tm(1) * 3600) - tm(2) * 60);%seconds

timeStr = num2str(tm(1));
if (tm(1) < 10)
    timeStr = ['0', timeStr];
end
timeStr = [timeStr, ':'];
if (tm(2) < 10)
    timeStr = [timeStr, '0'];
end
timeStr = [timeStr, num2str(tm(2))];
timeStr = [timeStr, ':'];
if (tm(3) < 10)
    timeStr = [timeStr, '0'];
end
timeStr = [timeStr, num2str(tm(3))];



function CreatDynamicFigs(hObject, eventdata, handles)
%creat figures after analysis
expAn = handles.ExpressAnalysis;%main form handle

hourEdtH = getappdata(expAn, 'hourEdtH');%handle of hour-edit
minEdtH = getappdata(expAn, 'minEdtH');%handle of minute-edit
etalonEdtH = getappdata(expAn, 'etalonEdtH');%handle of etalon-record-edit
nulTimeFig = getappdata(expAn, 'nulTimeFig');%handle of dialog window
nulTm = [str2double(get(hourEdtH, 'String')), str2double(get(minEdtH, 'String')), str2double(get(etalonEdtH, 'String'))];
rmappdata(expAn, 'hourEdtH');%remove handle
rmappdata(expAn, 'minEdtH');%remove handle
rmappdata(expAn, 'etalonEdtH');%remove handle
close(nulTimeFig)%close dialog

if ((numel(nulTm) < 3) || any(isnan(nulTm))) %no null time asigned
    return;
end
t0 = nulTm(1) * 60 + nulTm(2);%injection time (minutes from day begin)

set(handles.ExportData, 'UserData', nulTm);%null time (injection)
fBase = getappdata(expAn, 'fBase');%get list of recordations
rList = get(handles.recList, 'Data');%list of wanted files
rList = vertcat(rList{:, 2});%flags of accepted files
rCh = getappdata(expAn, 'rCh');%principal channel
prgMult = (-1) * str2double(get(handles.Porog, 'String'));%threshold multiplier
load(fBase(nulTm(3)).matFile, 'lfpVar')
etPrgs = lfpVar;%etalon thresholds

k = find(fBase(1).srcFile(1:(end - 2)) == '\', 1, 'last');
z = find(fBase(1).srcFile(1:(k - 2)) == '\', 1, 'last');
pf = [fBase(1).srcFile(1:k), 'graph\'];%pathname for figure saving
flNm = fBase(1).srcFile((z + 1):(k - 1));%directory name

% %% 1 plot SEP amplitude
% figure, title('SEP amplitude'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
%     if any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) %variable exit
%         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak')
%         if ismember('zp', varInFl)
%             load(fBase(w).matFile, 'zp')
%             zavp = zp;
%         else
%             load(fBase(w).matFile, 'zavp')
%         end
%         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
%             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'hd')
%             if ~isempty(sepOnsetPeak(rCh, sw).r)
%                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
%             end
%         elseif ismember('sepOnsetPeak', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak', 'hd')
%             sepOnsetPeak2 = sepOnsetPeak;
%         elseif ismember('sepOnsetPeak2', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak2', 'hd')
%         end
%         %sepOnsetPeak = CorrectSepOP(expAn, sepOnsetPeak, hd, zavp);%correct sepOnsetPeak
%         segms = GetSegments(hd, zavp);%segments of evoked activity
%         sgmNum = size(segms, 1);%number of segments
%         tm = zeros(sgmNum, 1);%time vector
%         bb = zeros(sgmNum, 1);%vector of parameter values
%         for t = 1:sgmNum %run over all segments
%             sw = segms(t, 3);%number of sweep where syncro-event was detected
%             g = segms(t, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%             %tm(t) = ((hd.inTTL_timestamps.t(t, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
%             tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%             bb(t) = sepOnsetPeak2(rCh, sw).r(g, 3);%value of parameters (SEP amplitude)
%         end
%         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%         tm = tm(isfinite(bb)); bb = abs(bb(isfinite(bb)));
%         plot(tm, bb, '.')
%         plot(mean(tm), mean(bb), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     end
% end
% plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(gca, 'YLabel'), 'String', 'amplitude, \muV')
% set(get(gca, 'XLabel'), 'String', 'time, min')
% set(gcf, 'FileName', [pf, 'SEP_amplitude.emf'], 'Name', [flNm, ' SEP_amplitude'], 'NumberTitle', 'off')
% 
% %% 2 plot SEP onset latency
% figure, title('SEP onset latency'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
%     if any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) %variable exit
%         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak')
%         if ismember('zp', varInFl)
%             load(fBase(w).matFile, 'zp')
%             zavp = zp;
%         else
%             load(fBase(w).matFile, 'zavp')
%         end
%         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
%             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'hd')
%             if ~isempty(sepOnsetPeak(rCh, sw).r)
%                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
%             end
%         elseif ismember('sepOnsetPeak', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak', 'hd')
%             sepOnsetPeak2 = sepOnsetPeak;
%         elseif ismember('sepOnsetPeak2', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak2', 'hd')
%         end
%         %sepOnsetPeak = CorrectSepOP(expAn, sepOnsetPeak, hd, zavp);%correct sepOnsetPeak
%         segms = GetSegments(hd, zavp);%segments of evoked activity
%         sgmNum = size(segms, 1);%number of segments
%         tm = zeros(sgmNum, 1);%time vector
%         bb = zeros(sgmNum, 1);%vector of parameter values
%         for t = 1:sgmNum %run over all segments
%             sw = segms(t, 3);%number of sweep where syncro-event was detected
%             g = segms(t, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%             %tm(t) = ((hd.inTTL_timestamps.t(t, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
%             tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%             bb(t) = sepOnsetPeak2(rCh, sw).r(g, 1);%value of parameters (SEP latency)
%         end
%         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%         tm = tm(isfinite(bb)); bb = bb(isfinite(bb));
%         plot(tm, bb, '.')
%         plot(mean(tm), mean(bb), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     end
% end
% plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(gca, 'YLabel'), 'String', 'latency, ms')
% set(get(gca, 'XLabel'), 'String', 'time, min')
% set(gcf, 'FileName', [pf, 'SEP_onset_latency.emf'], 'Name', [flNm ' SEP_onset_latency'], 'NumberTitle', 'off')
% 
% %% 3 plot MUA during SEP
% figure, title('MUA during SEP'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
%     if any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) %variable exit
%         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak', 'spks')
%         if ismember('zp', varInFl)
%             load(fBase(w).matFile, 'zp')
%             zavp = zp;
%         else
%             load(fBase(w).matFile, 'zavp')
%         end
%         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
%             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'hd', 'spks')
%             if ~isempty(sepOnsetPeak(rCh, sw).r)
%                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
%             end
%         elseif ismember('sepOnsetPeak', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak', 'hd', 'spks')
%             sepOnsetPeak2 = sepOnsetPeak;
%         elseif ismember('sepOnsetPeak2', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak2', 'hd', 'spks')
%         end
%         %sepOnsetPeak = CorrectSepOP(expAn, sepOnsetPeak, hd, zavp);%correct sepOnsetPeak
%         
%         segms = GetSegments(hd, zavp);%segments of evoked activity
%         sgmNum = size(segms, 1);%number of segments
%         tm = zeros(sgmNum, 1);%time vector
%         bb = zeros(sgmNum, 1);%vector of parameter values
%         for t = 1:sgmNum %run over all segments
%             sw = segms(t, 3);%number of sweep where syncro-event was detected
%             g = segms(t, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%             %tm(t) = ((hd.inTTL_timestamps.t(t, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
%             tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%             bb(t) = sepOnsetPeak2(rCh, sw).r(g, 4);%MUA during SEP
%         end
%         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%         tm = tm(isfinite(bb)); bb = bb(isfinite(bb));
%         plot(tm, bb, '.')
%         plot(mean(tm), mean(bb), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     end
% end
% plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(gca, 'YLabel'), 'String', 'MUA, units/20ms')
% set(get(gca, 'XLabel'), 'String', 'time, min')
% set(gcf, 'FileName', [pf, 'MUA_during_SEP.emf'], 'Name', [flNm, ' MUA_during_SEP'], 'NumberTitle', 'off')
% 
% %% 4 plot MUA after SEP
% figure, title('MUA after SEP'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
%     if any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) %variable exit
%         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak', 'spks')
%         if ismember('zp', varInFl)
%             load(fBase(w).matFile, 'zp')
%             zavp = zp;
%         else
%             load(fBase(w).matFile, 'zavp')
%         end
%         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
%             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'hd', 'spks')
%             if ~isempty(sepOnsetPeak(rCh, sw).r)
%                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
%             end
%         elseif ismember('sepOnsetPeak', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak', 'hd', 'spks')
%             sepOnsetPeak2 = sepOnsetPeak;
%         elseif ismember('sepOnsetPeak2', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak2', 'hd', 'spks')
%         end
%         %sepOnsetPeak = CorrectSepOP(expAn, sepOnsetPeak, hd, zavp);%correct sepOnsetPeak
%         
%         segms = GetSegments(hd, zavp);%segments of evoked activity
%         sgmNum = size(segms, 1);%number of segments
%         tm = zeros(sgmNum, 1);%time vector
%         bb = zeros(sgmNum, 1);%vector of parameter values
%         for t = 1:sgmNum %run over all segments
%             sw = segms(t, 3);%number of sweep where syncro-event was detected
%             g = segms(t, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%             %tm(t) = ((hd.inTTL_timestamps.t(t, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
%             tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%             bb(t) = sepOnsetPeak2(rCh, sw).r(g, 5);%MUA after SEP
%         end
%         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%         tm = tm(isfinite(bb)); bb = bb(isfinite(bb));
%         plot(tm, bb, '.')
%         plot(mean(tm), mean(bb), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     end
% end
% plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(gca, 'YLabel'), 'String', 'MUA, units/500ms')
% set(get(gca, 'XLabel'), 'String', 'time, min')
% set(gcf, 'FileName', [pf, 'MUA_after_SEP.emf'], 'Name', [flNm, ' MUA_after_SEP'], 'NumberTitle', 'off')
% 
% %% 5 plot Spontaneous MUA
% figure, title('Spontan MUA'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     load(fBase(w).matFile, 'zavp', 'hd', 'spks')
%     for sw = 1:hd.lActualEpisodes %run over sweeps
%         jj = spks(rCh, sw).ampl <= (prgMult * min(etPrgs(rCh, :)));
%         spks(rCh, sw).tStamp = spks(rCh, sw).tStamp(jj);
%         %spks(rCh, sw).ampl = spks(rCh, sw).ampl(jj);
%         %spks(rCh, sw).shape = spks(rCh, sw).shape(:, jj);
%     end
%     
%     if isempty(zavp.realStim.r) %~any(ttl > 0) %spontaneous
%         tm = 1:2:floor(diff(hd.recTime) / 60);%time vector
%         if isempty(tm) %too short recordation
%             tm = 1:2;%artificial time vector
%         end
%         spkDens = zeros(length(tm), 1);
%         sw = 1;
%         for z = 1:length(tm)
%             jj = ((tm(z) * 60) + [-150, 150]) * 1e3;
%             spkDens(z) = sum((spks(rCh, sw).tStamp > jj(1)) & (spks(rCh, sw).tStamp < jj(2)));
%         end
%         spkDens = spkDens / (diff(tm(1:2)) * 60);
%         tm = tm + ((hd.recTime(1) / 60) - t0);
%         plot(tm, spkDens, '.')
%         bins = [bins, repmat(tm(end), 2, 1)];
%         plot(mean(tm), mean(spkDens), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     else %evoked
%         segms = GetSegments(hd, zavp);%segments of evoked activity
%         sgmNum = size(segms, 1);%number of segments
%         tm = zeros(sgmNum, 1);%time vector
%         bb = zeros(sgmNum, 1);%vector of parameter values
%         for t = 1:sgmNum %run over all segments
%             sw = segms(t, 3);%number of sweep where syncro-event was detected
%             g = segms(t, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%             %tm(t) = ((hd.inTTL_timestamps.t(t, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
%             tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%             if (t > 1)
%                 spntMUA = sum((spks(rCh, sw).tStamp >= (segms(t - 1, 1) + 1000)) & (spks(rCh, sw).tStamp <= (segms(t, 1) - 0)));%MUA 1s after stim
%                 spntMUA = spntMUA * 1e3 / ((segms(t, 1) - 0) - (segms(t - 1, 1) + 1000));%unit/s
%             else %(z == 1) %first stimulus
%                 spntMUA = sum((spks(rCh, sw).tStamp >= 0) & (spks(rCh, sw).tStamp <= (segms(t, 1) - 0)));
%                 spntMUA = spntMUA * 1e3 / (segms(t, 1) - 0);%unit/s
%             end
%             bb(t) = spntMUA;%spontaneouse spikes per secound
%         end
%         plot(tm,  bb, '.')
%         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%         plot(mean(tm), mean(bb), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     end
% end
% plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(gca, 'YLabel'), 'String', 'MUA frequency, 1/s')
% set(get(gca, 'XLabel'), 'String', 'time, min')
% set(gcf, 'FileName', [pf, 'spontan_MUA.emf'], 'Name', [flNm, ' spontan_MUA'], 'NumberTitle', 'off')
% 
% %% 6-7 plot Evoked gamma- and alpha-beta range oscillations power
% fG = figure; axG = gca; title('Evoked \gamma-oscillations power'), hold on
% fAB = figure; axAB = gca; title('Evoked \alpha-\beta-oscillations power'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
%     if any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) %variable exit
%         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak', 'lfp')
%         if ismember('zp', varInFl)
%             load(fBase(w).matFile, 'zp')
%             zavp = zp;
%         else
%             load(fBase(w).matFile, 'zavp')
%         end
%         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
%             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'hd', 'lfp')
%             if ~isempty(sepOnsetPeak(rCh, sw).r)
%                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
%             end
%         elseif ismember('sepOnsetPeak', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak', 'hd', 'lfp')
%             sepOnsetPeak2 = sepOnsetPeak;
%         elseif ismember('sepOnsetPeak2', varInFl)
%             load(fBase(w).matFile, 'sepOnsetPeak2', 'hd', 'lfp')
%         end
%         lfp = double(lfp);
%         segms = GetSegments(hd, zavp);%segments of evoked activity
% 
%         segmEdge = [30, 530];%left and right shifts from synchro-point (ms)
%         whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to stimulus moments
%         whtLFP = squeeze(whtLFP);
% 
%         zavp.prm.Fs = 1 / (zavp.siS * zavp.rarStep);%sampling frequency for rare data
%         zavp.prm.fpass = [5, 90];
%         zavp.prm.tapers = [3, 5];%tapers
%         zavp.prm.pad = 2;%padding factors
%         zavp.prm.trialave = 0;%average over trials/channels when 1
%         [pwrSpct, frq] = mtspectrumc(whtLFP, zavp.prm);%mean power spectrum (mean by sweeps)
% 
%         sgmNum = size(segms, 1);%number of segments
%         tm = zeros(sgmNum, 1);%time vector
%         bb = zeros(sgmNum, 2);%vector of parameter values
%         for t = 1:sgmNum %run over all segments
%             %sw = segms(t, 3);%number of sweep where syncro-event was detected
%             %g = segms(t, 4);%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%             [~, p1] = min(abs(frq - 8));%low frequency
%             [~, p2] = min(abs(frq - 30));%high frequency
%             bb(t, 1) = sum(pwrSpct(p1:p2, t)) / (p2 - p1);%alpha-beta range
%             [~, p1] = min(abs(frq - 30));%low frequency
%             [~, p2] = min(abs(frq - 80));%high frequency
%             bb(t, 2) = sum(pwrSpct(p1:p2, t)) / (p2 - p1);%gamma range
%             %tm(t) = ((hd.inTTL_timestamps.t(t, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
%             tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%         end
%         plot(axAB, tm, bb(:, 1), 'b.')%alpha-beta range
%         plot(axG, tm, bb(:, 2), 'b.')%gamma range
%         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%         plot(axAB, mean(tm), mean(bb(:, 1)), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%         plot(axG, mean(tm), mean(bb(:, 2)), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     end
% end
% plot(axG, bins, repmat(ylim(axG)', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(axG, 'YLabel'), 'String', '\gamma-range power, \muV^2/Hz')
% set(get(axG, 'XLabel'), 'String', 'time, min')
% set(fG, 'FileName', [pf, 'Evk_gamma_oscillation_power.emf'], 'Name', [flNm, ' Evk_gamma_oscillation_power'], 'NumberTitle', 'off')
% 
% plot(axAB, bins, repmat(ylim(axAB)', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(axAB, 'YLabel'), 'String', '\alpha-\beta-range power, \muV^2/Hz')
% set(get(axAB, 'XLabel'), 'String', 'time, min')
% set(fAB, 'FileName', [pf, 'Evk_alpha-beta_oscillation_power.emf'], 'Name', [flNm, ' Evk_alpha-beta_oscillation_power'], 'NumberTitle', 'off')
% 
% 
% %% 8-9 plot Spontaneous gamma- and alpha-beta range oscillations power
% fG = figure; axG = gca; title('Spontan \gamma-oscillations power'), hold on
% fAB = figure; axAB = gca; title('Spontan \alpha-\beta-oscillations power'), hold on
% bins = [];%recordations begins
% for w = 1:length(fBase) %run over recordations
%     load(fBase(w).matFile, 'zavp', 'hd', 'lfp')
%     lfp = double(lfp);
%     if isempty(zavp.realStim.r) %~any(ttl > 0) %spontaneous
%         segms = zeros(floor(size(lfp, 1) / 1e3), 4);%preallocation of memory for stimuli moments
%         sgmNum = size(segms, 1);%number of segments
%         for z = sgmNum:-1:1 %
%             segms(z, 1) = size(lfp, 1) - ((sgmNum - z) * 1e3);%position of synchro-event (ms (from 0))
%             segms(z, 2) = zavp.stimCh;%number of channel where syncro-event was detected
%             segms(z, 3) = sw;%number of sweep where syncro-event was detected
%             segms(z, 4) = z;%number of stimulus (in matrix zavp.realStim) or number of trough in brst matrix
%         end
%         segmEdge = [-1e3, 0];%left and right shifts from synchro-point (ms)
%         whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to stimulus moments
%         whtLFP = squeeze(whtLFP);
%     else %evoked
%         segms = GetSegments(hd, zavp);%segments of evoked activity
%         segmEdge = [-1e3, 0];%left and right shifts from synchro-point (ms)
%         whtLFP = ZavSynchLFP(zavp, hd, segms, segmEdge, lfp, rCh, false);%lfp phased with respect to stimulus moments
%         whtLFP = squeeze(whtLFP);
%     end
%     zavp.prm.Fs = 1 / (zavp.siS * zavp.rarStep);%sampling frequency for rare data
%     zavp.prm.fpass = [5, 90];
%     zavp.prm.tapers = [3, 5];%tapers
%     zavp.prm.pad = 2;%padding factors
%     zavp.prm.trialave = 0;%average over trials/channels when 1
%     [pwrSpct, frq] = mtspectrumc(whtLFP, zavp.prm);%mean power spectrum (mean by sweeps)
%     
%     sgmNum = size(segms, 1);%number of segments
%     tm = zeros(sgmNum, 1);%time vector
%     bb = zeros(sgmNum, 2);%vector of parameter values
%     for t = 1:sgmNum %run over all segments
%         [~, p1] = min(abs(frq - 8));%low frequency
%         [~, p2] = min(abs(frq - 30));%high frequency
%         bb(t, 1) = sum(pwrSpct(p1:p2, t)) / (p2 - p1);%alpha-beta range
%         [~, p1] = min(abs(frq - 30));%low frequency
%         [~, p2] = min(abs(frq - 80));%high frequency
%         bb(t, 2) = sum(pwrSpct(p1:p2, t)) / (p2 - p1);%gamma range
%         %tm(t) = (segms(t, 1) * 1e-3 / 60) + ((hd.recTime(1) / 60) - t0);
%         tm(t) = hd.recTime(1) - t0 + segms(t, 1);%time after appl.
%     end
%     plot(axG, tm, bb(:, 2), 'b.')%gamma range
%     plot(axAB, tm, bb(:, 1), 'b.')%alpha-beta range
%     bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
%     plot(axG, mean(tm), mean(bb(:, 2)), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
%     plot(axAB, mean(tm), mean(bb(:, 1)), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
% end
% plot(axG, bins, repmat(ylim(axG)', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(axG, 'YLabel'), 'String', '\gamma-range power, \muV^2/Hz')
% set(get(axG, 'XLabel'), 'String', 'time, min')
% set(fG, 'FileName', [pf, 'Spnt_gamma_oscillation_power.emf'], 'Name', [flNm, ' Spnt_gamma_oscillation_power'], 'NumberTitle', 'off')
% 
% plot(axAB, bins, repmat(ylim(axAB)', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% set(get(axAB, 'YLabel'), 'String', '\alpha-\beta-range power, \muV^2/Hz')
% set(get(axAB, 'XLabel'), 'String', 'time, min')
% set(fAB, 'FileName', [pf, 'Spnt_alpha-beta_oscillation_power.emf'], 'Name', [flNm, ' Spnt_alpha-beta_oscillation_power'], 'NumberTitle', 'off')
% 
% % %% 7 plot first-second trough
% % figure, title('first-second trough'), hold on
% % bins = [];%recordations begins
% % for w = 1:length(fBase) %run over recordations
% %     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
% %     if (any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) && ismember('scndTrg', varInFl)) %variable exit
% %         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak', 'scndTrg')
% %         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
% %             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'zavp', 'hd')
% %             if ~isempty(sepOnsetPeak(rCh, sw).r)
% %                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
% %             end
% %         elseif ismember('sepOnsetPeak', varInFl)
% %             load(fBase(w).matFile, 'sepOnsetPeak', 'zavp', 'hd')
% %             sepOnsetPeak2 = sepOnsetPeak;
% %         elseif ismember('sepOnsetPeak2', varInFl)
% %             load(fBase(w).matFile, 'sepOnsetPeak2', 'zavp', 'hd')
% %         end
% %         ii = all(isfinite(sepOnsetPeak(rCh, sw).r(:, 1:3)), 2) & all(isfinite(scndTrg(1).r(:, 1:3)), 2);
% %         %tm = ((hd.inTTL_timestamps.t(:, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
% %         tm = hd.recTime(1) - t0 + segms(:, 1);%time after appl.
% %         plot(tm(ii), scndTrg(1).r(ii, 3) ./ sepOnsetPeak(rCh, sw).r(ii, 3), '.')
% %         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
% %         plot(mean(tm(ii)), mean(abs(scndTrg(1).r(ii, 3) ./ sepOnsetPeak(rCh, sw).r(ii, 3))), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
% %     end
% % end
% % plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% % set(get(gca, 'YLabel'), 'String', 'scnd/first ratio')
% % set(gcf, 'FileName', [pf, 'first-second_trough.emf'], 'Name', [flNm, ' first-second_trough'], 'NumberTitle', 'off')
% % 
% % %% 8 second peak latency
% % figure, title('second peak latency'), hold on
% % bins = [];%recordations begins
% % for w = 1:length(fBase) %run over recordations
% %     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
% %     if (any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) && ismember('scndTrg', varInFl)) %variable exit
% %         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak', 'scndTrg')
% %         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
% %             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'zavp', 'hd')
% %             if ~isempty(sepOnsetPeak(rCh, sw).r)
% %                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
% %             end
% %         elseif ismember('sepOnsetPeak', varInFl)
% %             load(fBase(w).matFile, 'sepOnsetPeak', 'zavp', 'hd')
% %             sepOnsetPeak2 = sepOnsetPeak;
% %         elseif ismember('sepOnsetPeak2', varInFl)
% %             load(fBase(w).matFile, 'sepOnsetPeak2', 'zavp', 'hd')
% %         end
% %         ii = all(isfinite(sepOnsetPeak(rCh, sw).r(:, 1:3)), 2) & all(isfinite(scndTrg(1).r(:, 1:3)), 2);
% %         %tm = ((hd.inTTL_timestamps.t(:, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
% %         tm = hd.recTime(1) - t0 + segms(:, 1);%time after appl.
% %         plot(tm(ii), scndTrg(1).r(ii, 2), '.')
% %         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
% %         plot(mean(tm(ii)), mean(scndTrg(1).r(ii, 2)), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
% %     end
% % end
% % plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% % set(get(gca, 'YLabel'), 'String', 'scnd peak latency, ms')
% % set(gcf, 'FileName', [pf, 'second_peak_latency.emf'], 'Name', [flNm, 'second_peak_latency'], 'NumberTitle', 'off')
% % 
% % %% 9 plot MUA2t/MUA1t
% % figure, title('MUA2t/MUA1t'), hold on
% % bins = [];%recordations begins
% % for w = 1:length(fBase) %run over recordations
% %     varInFl = who(matfile(fBase(w).matFile));%variables in requested file
% %     if (any(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl)) && ismember('scndTrg', varInFl)) %variable exit
% %         %load(fBase(w).matFile, 'zavp', 'hd', 'sepOnsetPeak', 'scndTrg')
% %         if all(ismember({'sepOnsetPeak', 'sepOnsetPeak2'}, varInFl))
% %             load(fBase(w).matFile, 'sepOnsetPeak', 'sepOnsetPeak2', 'zavp', 'hd')
% %             if ~isempty(sepOnsetPeak(rCh, sw).r)
% %                 sepOnsetPeak2(rCh, sw).r(:, 1:5) = sepOnsetPeak(rCh, sw).r(:, 1:5);%replase by manually corrected data
% %             end
% %         elseif ismember('sepOnsetPeak', varInFl)
% %             load(fBase(w).matFile, 'sepOnsetPeak', 'zavp', 'hd')
% %             sepOnsetPeak2 = sepOnsetPeak;
% %         elseif ismember('sepOnsetPeak2', varInFl)
% %             load(fBase(w).matFile, 'sepOnsetPeak2', 'zavp', 'hd')
% %         end
% %         ii = all(isfinite(sepOnsetPeak(rCh, sw).r(:, 1:3)), 2) & all(isfinite(scndTrg(1).r(:, 1:3)), 2);
% %         %tm = ((hd.inTTL_timestamps.t(:, 1) * 1e-6) / 60) + ((hd.recTime(1) / 60) - t0);
% %         tm = hd.recTime(1) - t0 + segms(:, 1);%time after appl.
% %         plot(tm(ii), scndTrg(1).r(ii, 4) ./ sepOnsetPeak(rCh, sw).r(ii, 4), '.')
% %         bins = [bins, repmat((hd.recTime(1) / 60) - t0, 2, 1)];%add time of begin of current recordation
% %         jj = scndTrg(1).r(ii, 4) ./ sepOnsetPeak(rCh, sw).r(ii, 4);
% %         plot(mean(tm(ii)), mean(jj(isfinite(jj))), 'or', 'MarkerFaceColor', [1, 0, 0], 'MarkerEdgeColor', [0, 0, 0], 'MarkerSize', 8)
% %     end
% % end
% % plot(bins, repmat(ylim', 1, size(bins, 2)), '--k', 'LineWidth', 2)
% % set(get(gca, 'YLabel'), 'String', 'MUA2t / MUA1t')
% % set(gcf, 'FileName', [pf, 'MUA2t-MUA1t.emf'], 'Name', [flNm, ' MUA2t-MUA1t'], 'NumberTitle', 'off')


% --- Executes when selected cell(s) is changed in recList.
function recList_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to recList (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)

ii = unique(eventdata.Indices(:, 1));%indices of selected rows
if ((length(ii) > 1) && all(diff(ii) == 1) && all(eventdata.Indices(:, 2) == 2)) %change channels selection
    recListData = get(handles.recList, 'Data');%read current table content
    jj = horzcat(recListData{ii(:, 1), 2});%selection falgs
    trg = ~(sum(jj) >= sum(~jj));%target value (prevailing selection)
    for z = ii(:, 1)' %run over raws
        recListData{z, 2} = trg;%change value of selected cells (set prevailing value)
    end
    set(handles.recList, 'Data', recListData)%set table data
end
