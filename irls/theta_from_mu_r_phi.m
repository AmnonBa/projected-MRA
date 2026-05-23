function theta = theta_from_mu_r_phi(T1hat, r, phi, p)
    q = (p - 1)/2;

    thetaHat0 = sqrt(p)/2 * mean(T1hat);

    thetaHat = zeros(p,1);
    thetaHat(1) = thetaHat0;

    for k = 1:q
        thetaHat(k+1) = r(k) * exp(1i*phi(k));
        thetaHat(p-k+1) = conj(thetaHat(k+1));
    end

    theta = real(ifft(thetaHat * sqrt(p)));
end

