function [sBest, kBest] = nearest_sheet(x, beta)
    bestVal = inf;
    sBest = 1;
    kBest = 0;

    for s = [-1, 1]
        k = round((x - s*beta) / (2*pi));
        val = abs(x - s*beta - 2*pi*k);
        if val < bestVal
            bestVal = val;
            sBest = s;
            kBest = k;
        end
    end
end
