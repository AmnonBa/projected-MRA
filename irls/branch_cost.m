function c = branch_cost(eps, beta, dHat, dStar, q)
    residuals = [];

    for j = 2:q-2
        pred = cos(-eps(1)*beta(1) + eps(j)*beta(j) + eps(j+1)*beta(j+1));
        residuals = [residuals; pred - dHat(j-1)]; 
    end

    predStar = cos(-eps(1)*beta(1) + eps(q-1)*beta(q-1) + eps(q)*beta(q));
    residuals = [residuals; predStar - dStar];

    c = sum(residuals.^2);
end
