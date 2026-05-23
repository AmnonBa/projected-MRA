function [betaChain, dHat, dStar, valid] = extract_chain_quantities(M3hat, rHat, q, clipCosines)
    valid = true;
    betaChain = zeros(q,1);

    c = zeros(q,1);

    c(1) = M3hat(1,1,2) / (2*rHat(1)^2*rHat(2));

    for j = 2:q-1
        c(j) = M3hat(1,j,j+1) / (2*rHat(1)*rHat(j)*rHat(j+1));
    end

    c(q) = M3hat(1,q,q) / (2*rHat(1)*rHat(q)^2);

    if clipCosines
        c = max(min(real(c),1),-1);
    elseif any(abs(c) > 1) || any(~isfinite(c))
        valid = false;
        dHat = [];
        dStar = [];
        return;
    end

    betaChain = acos(c);

    if q >= 4
        dHat = zeros(q-3,1);
        for j = 2:q-2
            dHat(j-1) = M3hat(2,j,j+2) / (2*rHat(2)*rHat(j)*rHat(j+2));
        end
    else
        dHat = [];
    end

    dStar = M3hat(2,q-1,q) / (2*rHat(2)*rHat(q-1)*rHat(q));

    if clipCosines
        dHat = max(min(real(dHat),1),-1);
        dStar = max(min(real(dStar),1),-1);
    end
end
