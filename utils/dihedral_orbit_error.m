function err = dihedral_orbit_error(thetaEst, thetaTrue)
    p = numel(thetaTrue);

    thetaEst = thetaEst(:);
    thetaTrue = thetaTrue(:);

    thetaRef = reflect_signal(thetaTrue);

    bestErr = inf;
    for ell = 0:p-1
        cand1 = apply_cyclic_shift(thetaTrue, ell);
        cand2 = apply_cyclic_shift(thetaRef, ell);

        bestErr = min(bestErr, norm(thetaEst - cand1));
        bestErr = min(bestErr, norm(thetaEst - cand2));
    end

    err = bestErr;
end
