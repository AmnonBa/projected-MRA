function [M2hat, M3hat] = transfer_to_cosine_moments(T1hat, T2hat, T3hat, Ainv, p)
    q = numel(T1hat);
    mu = T1hat;

    T2center = T2hat - mu*mu.';

    T3center = zeros(q,q,q);
    for a = 1:q
        for b = 1:q
            for c = 1:q
                T3center(a,b,c) = ...
                    T3hat(a,b,c) ...
                    - mu(a)*T2hat(b,c) ...
                    - mu(b)*T2hat(a,c) ...
                    - mu(c)*T2hat(a,b) ...
                    + 2*mu(a)*mu(b)*mu(c);
            end
        end
    end

    M2hat = p * Ainv * T2center * Ainv.';

    M3hat = zeros(q,q,q);
    scale = p^(3/2);

    for k1 = 1:q
        for k2 = 1:q
            for k3 = 1:q
                val = 0;
                for j1 = 1:q
                    for j2 = 1:q
                        for j3 = 1:q
                            val = val + Ainv(k1,j1) * Ainv(k2,j2) * Ainv(k3,j3) ...
                                * T3center(j1,j2,j3);
                        end
                    end
                end
                M3hat(k1,k2,k3) = scale * val;
            end
        end
    end
end
