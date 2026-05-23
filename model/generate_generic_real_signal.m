function theta = generate_generic_real_signal(p, minFourierMag)
    q = (p - 1)/2;
    theta = randn(p,1);
    theta = theta / norm(theta);
    thetaHat = fft(theta) / sqrt(p);

    while min(abs(thetaHat(2:q+1))) < minFourierMag
        theta = randn(p,1);
        theta = theta / norm(theta);
        thetaHat = fft(theta) / sqrt(p);
    end
end
