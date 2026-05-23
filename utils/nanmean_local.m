function m = nanmean_local(x)
    x = x(:);
    x = x(~isnan(x));
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end
