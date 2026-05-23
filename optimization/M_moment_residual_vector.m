function r = M_moment_residual_vector(theta, T1hat, M2hat, M3hat)
    [T1, ~, ~] = projected_moments_from_theta(theta);
    [M2, M3] = cosine_moments_from_theta(theta);

    s1 = max(norm(T1hat(:)), 1e-8);
    s2 = max(norm(M2hat(:)), 1e-8);
    s3 = max(norm(M3hat(:)), 1e-8);

    r1 = (T1(:) - T1hat(:)) / s1;
    r2 = (M2(:) - M2hat(:)) / s2;
    r3 = (M3(:) - M3hat(:)) / s3;

    r = [r1; r2; r3];
end


