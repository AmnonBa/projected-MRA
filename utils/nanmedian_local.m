function m = nanmedian_local(x)
    x = x(:);
    x = x(~isnan(x));
    if isempty(x)
        m = NaN;
    else
        m = median(x);
    end
end
