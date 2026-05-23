function s = mean_success(E, tol)
    s = zeros(size(E,1),1);
    for i = 1:size(E,1)
        x = E(i,:);
        x = x(~isnan(x));
        if isempty(x)
            s(i) = NaN;
        else
            s(i) = mean(x < tol);
        end
    end
end
