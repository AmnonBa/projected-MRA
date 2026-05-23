function epsBest = viterbi_sign_recovery(beta, dHat, dStar, q, eps1)
    states = [-1, +1];

    if q == 3
        bestCost = inf;
        epsBest = [];
        for e2 = states
            for e3 = states
                eps = [eps1; e2; e3];
                c = branch_cost(eps, beta, dHat, dStar, q);
                if c < bestCost
                    bestCost = c;
                    epsBest = eps;
                end
            end
        end
        return;
    end

    dp = inf(q,2);
    prev = zeros(q,2);
    dp(2,:) = 0;

    for j = 2:q-2
        for sIdx = 1:2
            s = states(sIdx);
            for tIdx = 1:2
                t = states(tIdx);
                pred = cos(-eps1*beta(1) + s*beta(j) + t*beta(j+1));
                obs = dHat(j-1);
                localCost = (pred - obs)^2;
                newCost = dp(j,sIdx) + localCost;

                if newCost < dp(j+1,tIdx)
                    dp(j+1,tIdx) = newCost;
                    prev(j+1,tIdx) = sIdx;
                end
            end
        end
    end

    bestCost = inf;
    bestPrevIdx = 1;
    bestLastIdx = 1;

    for sIdx = 1:2
        s = states(sIdx);
        for tIdx = 1:2
            t = states(tIdx);
            pred = cos(-eps1*beta(1) + s*beta(q-1) + t*beta(q));
            totalCost = dp(q-1,sIdx) + (pred - dStar)^2;

            if totalCost < bestCost
                bestCost = totalCost;
                bestPrevIdx = sIdx;
                bestLastIdx = tIdx;
            end
        end
    end

    epsBest = zeros(q,1);
    epsBest(1) = eps1;
    epsBest(q) = states(bestLastIdx);
    epsBest(q-1) = states(bestPrevIdx);

    curIdx = bestPrevIdx;
    for j = q-1:-1:3
        prevIdx = prev(j,curIdx);
        epsBest(j-1) = states(prevIdx);
        curIdx = prevIdx;
    end
end
