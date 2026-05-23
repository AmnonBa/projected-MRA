function [theta, llTrace, info] = em_projected_mra(Y, theta0, sigmaNoise, maxIters, tol)
    [n,q] = size(Y);
    p = 2*q + 1;

    theta = theta0(:);
    llTrace = zeros(maxIters,1);

    Xall = zeros(p,q);

    info.iterations = 0;
    info.converged = false;
    info.finalRelChange = NaN;

    for iter = 1:maxIters
        for ell = 0:p-1
            Xall(ell+1,:) = projection_pi(apply_cyclic_shift(theta, ell)).';
        end

        logW = zeros(n,p);
        for ell = 1:p
            diff = Y - Xall(ell,:);
            sqdist = sum(diff.^2,2);
            logW(:,ell) = -0.5 * sqdist / sigmaNoise^2;
        end

        [logNorm, W] = normalize_log_weights(logW);
        llTrace(iter) = sum(logNorm) - n*log(p) - n*q*log(sigmaNoise*sqrt(2*pi));

        H = zeros(p,p);
        b = zeros(p,1);

        for ell = 0:p-1
            Pmat = projection_shift_matrix(p, ell);
            w = W(:,ell+1);
            wsum = sum(w);

            H = H + wsum * (Pmat.' * Pmat);
            b = b + Pmat.' * (Y.' * w);
        end

        lambda = 1e-10;
        thetaNew = (H + lambda*eye(p)) \ b;

        relChange = norm(thetaNew - theta) / max(1,norm(theta));
        theta = thetaNew;

        info.iterations = iter;
        info.finalRelChange = relChange;

        if relChange < tol
            info.converged = true;
            llTrace = llTrace(1:iter);
            return;
        end
    end

    llTrace = llTrace(1:maxIters);
end
