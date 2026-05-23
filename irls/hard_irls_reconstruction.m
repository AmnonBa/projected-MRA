function [thetaIRLS, info] = hard_irls_reconstruction(T1hat, M2hat, M3hat, p, irls)
    q = (p - 1)/2;

    info.iterations = 0;
    info.converged = false;
    info.finalRelChange = NaN;

    rHat = zeros(q,1);
    for k = 1:q
        rHat(k) = sqrt(max(real(M2hat(k,k))/2, 1e-12));
    end

    [betaChain, dHat, dStar, valid] = extract_chain_quantities(M3hat, rHat, q, irls.clipCosines);

    if ~valid
        thetaIRLS = zeros(p,1);
        return;
    end

    epsInit = viterbi_sign_recovery(betaChain, dHat, dStar, q, +1);
    phi0 = phases_from_branch(betaChain, epsInit, p);

    eqs = build_all_phase_equations(M3hat, rHat, p, irls.clipCosines);

    if isempty(eqs.B)
        thetaIRLS = theta_from_mu_r_phi(T1hat, rHat, phi0, p);
        return;
    end

    phi = phi0;

    B = eqs.B;
    beta = eqs.beta;
    amp = eqs.amp;

    wBase = amp.^2;
    wBase = wBase / max(mean(wBase), 1e-12);

    for iter = 1:irls.numIters
        x = B * phi;
        b = zeros(size(beta));

        for rr = 1:numel(beta)
            [sBest, kBest] = nearest_sheet(x(rr), beta(rr));
            b(rr) = sBest * beta(rr) + 2*pi*kBest;
        end

        Wsqrt = sqrt(wBase(:));
        Wmat = spdiags(Wsqrt, 0, numel(Wsqrt), numel(Wsqrt));

        Bweighted = Wmat * B;
        bweighted = Wsqrt .* b;

        Baug = [Bweighted; sqrt(irls.gaugeWeight)*speye(q)];
        baug = [bweighted; sqrt(irls.gaugeWeight)*phi0];

        phiNew = Baug \ baug;

        relChange = norm(phiNew - phi) / max(1, norm(phi));
        phi = phiNew;

        info.iterations = iter;
        info.finalRelChange = relChange;

        if iter >= irls.minIters && relChange < irls.tol
            info.converged = true;
            break;
        end
    end

    thetaIRLS = theta_from_mu_r_phi(T1hat, rHat, phi, p);
end
