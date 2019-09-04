function [pntsL, pntsR] = ZavFindFlats(data)
%stSt = ZavFindFlats(data)
%looking for "flat" segments (uninterrupted segments with same value)
%
%INPUTS
%data - array of numbers
%
%OUTPUS
%pntsL - left points of flat segments
%pntsR - rigth points of flat segments

pntsL = find(abs(diff(data)) < 1e-10);%left points of flat segments
pntsR = pntsL(end:-1:1);%rigth points of flat segments

reps = find(diff(pntsL) == 1);%neighbor indices
while ~isempty(reps)
    pntsL(reps + 1) = [];%exclude repetitions
    reps = find(diff(pntsL) == 1);%find repetitions again
end

reps = find(diff(pntsR) == -1);%neighbor indices
while ~isempty(reps)
    pntsR(reps + 1) = [];%exclude repetitions
    reps = find(diff(pntsR) == -1);%find repetitions again
end
pntsR = sort(pntsR) + 1;%sort right points
