function [spk, ampl] = ZavFindSpikes(dataFlt, fs, thrshld, fltFlg, durSort)
%spk = ZavFindSpikes(dataFlt, fs, thrshld)
%looking for spikes on single trace (for a sweep on a channel)
%
%INPUTS
%dataFlt - single trace (sweep on a channel (single dimensional)
%fs - sampling frequency (Hz)
%thrshld - threshold defined by user (enter [], if not need)
%fltFlg - filtration flag (do filtration if true)
%durSort - do sorting by duration if true
%
%OUTPUTS
%spk - moments of detected spikes (samples)
%ampl - amplitude of spikes

if (fltFlg) %do filtration if true
    %highpass filtration
%     dataFlt = ZavFilter(dataFlt, fs, 'high', 400, 2);%filtration (mode 2 - cascade filtration)
    
    %bandpass filtration
    dataFlt = ZavFilter(dataFlt, fs, 'high', 300, 2);%highpass filtration
    if ((fs / 2.001) > 5000) %fast discretisation
        dataFlt = ZavFilter(dataFlt, fs, 'low', 5000, 2);%lowpass filtration
    end
end
if isempty(thrshld)
    thrshld = -3 * std(dataFlt, 0, 1);%threshold
end

% %%% method 1: domains of connectedness %%%
% spk = find(dataFlt < thrshld);%amplitude criterion
% reps = find(diff(spk) == 1);%repetitions on index
% while ~isempty(reps)
%     spk(reps + 1) = [];%exclude repetitions
%     reps = find(diff(spk) == 1);%find repetitions again
% end
% %%% end of method 1 %%%

%%% method 2: local minima %%%
spks1 = find(dataFlt < thrshld);%amplitude criterion
spks1((spks1 <= 1) | (spks1 >= numel(dataFlt))) = [];%delete edge element because of command (spks1(t1) -/+ 1) in next lines

spk = zeros(numel(spks1), 1);
z = 1;
interrupts = find(diff(spks1) > 1);%points of interruptions
t1 = 1;%left board of domain of connectedness
for t = 1:numel(interrupts) %run over all interruption points
    jj = (spks1(t1) - 1):(spks1(interrupts(t)) + 1);
    mins = ZavFindMins(dataFlt(jj));%multiple minima expected;%[~, mins] = min(dataFlt(jj));%for short segments
    spk(z:(z + numel(mins) - 1)) = mins + jj(1) - 1;
    z = z + numel(mins);
    t1 = interrupts(t) + 1;
end
if (t1 <= numel(spks1))
    mins = ZavFindMins(dataFlt((spks1(t1) - 1):(spks1(end) + 1)));
    spk(z:(z + numel(mins) - 1)) = mins + spks1(t1) - 2;
    z = z + numel(mins);
end
spk(z:end) = [];%delete excess elements
%%% end of method 2 %%%

%%% method 3: localminima %%%
%%% end of method 3 %%%

%sort spike by length (select spikes longer than 0.2 ms on threshold level)
if (durSort)
    %way 1: wide on threshold level
    t1 = round(0.2e-3 / (1 / fs));%length of segment (samples)
    t2 = round(0.3e-3 / (1 / fs));%minimal distance between spikes (samples) 0.3
    jj = zeros(1, 2);
    spk = spk((spk > (5 * t1)) & (spk <= (size(dataFlt, 1) - (5 * t1))));%remove spikes on edges
    mins = ones(numel(spk), 1);
    for t = 1:numel(spk)
        if (t > 1)
            if ((spk(t) - spk(t - 1)) < t2)
                mins(t) = 0;%remember number of wrong spike (too close to previous spike)
            end
        end
        if (mins(t) == 1)%not too close to previous spike
            spks1 = dataFlt(spk(t) + (0:-1:(-5 * t1)));%left side
            p1 = find(spks1 >= thrshld, 1, 'first');%first point under threshold
            if isempty(p1)
                jj(1) = t1;
            else
                jj(1) = interp1(dataFlt(spk(t) - ([p1 - 1, p1] - 1)), [p1 - 1, p1], thrshld) - 1;%count of points above threshold on left side
            end

            spks1 = dataFlt(spk(t) + (0:(5 * t1)));%right side
            p1 = find(spks1 >= thrshld, 1, 'first');%first point under threshold
            if isempty(p1)
                jj(2) = t1;
            else
                jj(2) = interp1(dataFlt(spk(t) + ([p1 - 1, p1] - 1)), [p1 - 1, p1], thrshld) - 1;%count of points above threshold on right side
            end

            if (sum(jj) < t1)%too short spike
                mins(t) = 0;%remember number of wrong spike (to short spike)
            end
        end
    end
    
%     %way 2: assessment of full wide
%     t1 = round(0.2e-3 / (1 / fs));%length of segment (samples)
%     t2 = round(0.3e-3 / (1 / fs));%minimal distance between spikes
%     jj = zeros(1, 2);
%     spk = spk((spk > (5 * t1)) & (spk <= (size(dataFlt, 1) - (5 * t1))));%remove spikes on edges
%     mins = ones(numel(spk), 1);
%     for t = 1:numel(spk)
%         if (t > 1)
%             if ((spk(t) - spk(t - 1)) < t2)
%                 mins(t) = 0;%remember number of wrong spike (too close to previous spike)
%             end
%         end
%         if (mins(t) == 1)%not too close to previous spike
%             spks1 = dataFlt(spk(t) + (0:-1:(-5 * t1)));%left side
%             p1 = find(diff(spks1) <= 0.1, 1, 'first');%left nearest maximum
%             if isempty(p1)
%                 jj(1) = t1;
%             else
%                 jj(1) = p1;
%             end
% 
%             spks1 = dataFlt(spk(t) + (0:(5 * t1)));%right side
%             p1 = find(diff(spks1) <= 0.1, 1, 'first');%right nearest maximum
%             if isempty(p1)
%                 jj(2) = t1;
%             else
%                 jj(2) = p1;
%             end
% 
%             if (sum(jj) < t1)%too short spike
%                 mins(t) = 0;%remember number of wrong spike (to short spike)
%             end
%         end
%     end
    
    spk = spk(mins == 1);%right spikes remain only
end
ampl = dataFlt(spk);
