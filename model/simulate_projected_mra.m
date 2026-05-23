function Y = simulate_projected_mra(theta, n, sigmaNoise)
    p = numel(theta);
    q = (p - 1)/2;
    Y = zeros(n,q);

    for i = 1:n
        ell = randi(p) - 1;
        shifted = apply_cyclic_shift(theta, ell);
        xProj = projection_pi(shifted);
        Y(i,:) = (xProj + sigmaNoise * randn(q,1)).';
    end
end
