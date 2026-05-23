function [logNorm, W] = normalize_log_weights(logW)
    m = max(logW, [], 2);
    expW = exp(logW - m);
    s = sum(expW,2);
    W = expW ./ s;
    logNorm = m + log(s);
end
