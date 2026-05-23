function r = T_moment_residual_vector(theta, T1hat, T2hat, T3hat)
    [T1, T2, T3] = projected_moments_from_theta(theta);

    s1 = max(norm(T1hat(:)), 1e-8);
    s2 = max(norm(T2hat(:)), 1e-8);
    s3 = max(norm(T3hat(:)), 1e-8);

    r1 = (T1(:) - T1hat(:)) / s1;
    r2 = (T2(:) - T2hat(:)) / s2;
    r3 = (T3(:) - T3hat(:)) / s3;

    r = [r1; r2; r3];
end
