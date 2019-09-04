function mins = ZavFindMins(y, frad, deep)
%mins = ZavFindMins(y, frad, deep)
%find local minima
%
%INPUTS
%y - signal (single dimensional)
%frad - frequency radius of minima
%deep - deep of minimum (amplitude units)
%
%OUTPUS
%mins - point of local minima

if (size(y, 1) == 1)
    y = y';
end

dY = diff(y);%first derivative

%find nulls of first derivative with positive second derivative (points of minimum)
if (~isempty(dY))
    %derivative signature is changed from negative to non-negative! (nulls of funtion dY)
    %dY_nuls = find((sign(dY) == 0) | ((sign([dY(2:end); dY(end)]) - sign(dY)) > 1));%primary nuls
    dY_nuls = find(diff(dY) >= 0);%primary nuls of dY
    mins = dY_nuls + 1;%correction of indices on shift due to derivation
    
    %direct check of minima
    prg = (max(y) - min(y)) / 1e10;%magnitude threshold for minima
    mins((mins < 2) | (mins >= numel(y))) = [];%delete minima on edges
    ii = true(size(mins));
    for t = 1:numel(mins)
        g = mins(t) + (-1:1);
        g = y(g);
        if ~((((g(1) - g(2)) > prg) && ((g(3) - g(2)) > prg)) || ...
             (((g(1) - g(2)) > prg) && (g(3) == g(2))) || ...
             (((g(3) - g(2)) > prg) && (g(1) == g(2))) ...
           )
            ii(t) = false;%delete false minimum
        end
    end
    mins = mins(ii);
    if (nargin > 1) %frequency radius of minima (frad) and deep of minimum was assigned
        %frad - frequency radius of minima
        sprWid = 7;%width of separating wall between trough
        mins = mins((mins > (frad * sprWid)) & (mins < (numel(y) - (frad * sprWid) + 1)));%non edge nulls
        shft = zeros(frad, sprWid);%array of shifts
        for t = 1:frad
            shft(t, :) = t:(t + sprWid - 1);%array of shifts
        end
        t = 1;%
        while (t < numel(mins))
            g = 1;
            leftF = true;%left flag of strong minimum
            rightF = true;%right flag of strong minimum
            while ((g <= frad) && (leftF || rightF))
                if (leftF && all((y(mins(t) - shft(g, :)) - y(mins(t))) >= deep))
                    leftF = false;%strong minimum (left side)
                end
                if (rightF && all((y(mins(t) + shft(g, :)) - y(mins(t))) >= deep))
                    rightF = false;%strong minimum (right side)
                end
                g = g + 1;
            end
            if (~leftF && ~rightF)%strong minimum
                t = t + 1;%accept strong minimum
            else
                mins(t) = [];%delete weak minimum
            end
        end
    end
else
    mins = [];%no minima
end
