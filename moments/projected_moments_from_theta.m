function [T1, T2, T3] = projected_moments_from_theta(theta)
    p = numel(theta);
    q = (p - 1)/2;

    X = zeros(q,p);
    for ell = 0:p-1
        X(:,ell+1) = projection_pi(apply_cyclic_shift(theta, ell));
    end

    T1 = mean(X,2);
    T2 = (X * X.') / p;

    T3 = zeros(q,q,q);
    for ell = 1:p
        x = X(:,ell);
        for a = 1:q
            for b = 1:q
                xab = x(a)*x(b);
                for c = 1:q
                    T3(a,b,c) = T3(a,b,c) + xab*x(c);
                end
            end
        end
    end
    T3 = T3 / p;
end
