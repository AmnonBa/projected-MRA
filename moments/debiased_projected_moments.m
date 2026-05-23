function [T1hat, T2hat, T3hat] = debiased_projected_moments(Y, sigmaNoise)
    [n,q] = size(Y);

    T1raw = mean(Y,1).';

    T2raw = zeros(q,q);
    T3raw = zeros(q,q,q);

    for i = 1:n
        y = Y(i,:).';

        T2raw = T2raw + y*y.';

        for a = 1:q
            for b = 1:q
                yab = y(a)*y(b);
                for c = 1:q
                    T3raw(a,b,c) = T3raw(a,b,c) + yab*y(c);
                end
            end
        end
    end

    T2raw = T2raw / n;
    T3raw = T3raw / n;

    T1hat = T1raw;
    T2hat = T2raw - sigmaNoise^2 * eye(q);

    T3hat = T3raw;
    for a = 1:q
        for b = 1:q
            for c = 1:q
                correction = 0;
                if a == b
                    correction = correction + sigmaNoise^2 * T1hat(c);
                end
                if a == c
                    correction = correction + sigmaNoise^2 * T1hat(b);
                end
                if b == c
                    correction = correction + sigmaNoise^2 * T1hat(a);
                end
                T3hat(a,b,c) = T3hat(a,b,c) - correction;
            end
        end
    end
end
