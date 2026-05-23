function [M2, M3] = cosine_moments_from_theta(theta)
    p = numel(theta);
    q = (p - 1)/2;

    thetaHat = fft(theta) / sqrt(p);

    C = zeros(p,q);
    for ell = 0:p-1
        for k = 1:q
            z = thetaHat(k+1) * exp(-2*pi*1i*k*ell/p);
            C(ell+1,k) = 2*real(z);
        end
    end

    M2 = (C.' * C) / p;

    M3 = zeros(q,q,q);
    for ell = 1:p
        cvec = C(ell,:).';
        for a = 1:q
            for b = 1:q
                cab = cvec(a)*cvec(b);
                for c = 1:q
                    M3(a,b,c) = M3(a,b,c) + cab*cvec(c);
                end
            end
        end
    end
    M3 = M3 / p;
end
