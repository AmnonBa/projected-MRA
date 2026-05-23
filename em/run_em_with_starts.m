function [thetaBest, bestLL, info] = run_em_with_starts(Y, sigmaNoise, thetaIRLS, em)
    [~, q] = size(Y);
    p = 2*q + 1;

    starts = {};

    if em.includeIRLSStart
        starts{end+1} = thetaIRLS(:); 
    end

    for s = 1:em.numRandomStarts
        z = randn(p,1);
        z = z / norm(z) * max(norm(thetaIRLS), 1);
        starts{end+1} = z; 
    end

    bestLL = -inf;
    thetaBest = starts{1};

    info.totalIterations = 0;
    info.bestIterations = NaN;
    info.numStarts = numel(starts);
    info.numSuccessfulStarts = 0;

    for s = 1:numel(starts)
        [thetaCand, llTrace, infoStart] = em_projected_mra(Y, starts{s}, sigmaNoise, em.numIters, em.tol);
        finalLL = llTrace(end);

        info.numSuccessfulStarts = info.numSuccessfulStarts + 1;
        info.totalIterations = info.totalIterations + infoStart.iterations;

        if finalLL > bestLL
            bestLL = finalLL;
            thetaBest = thetaCand;
            info.bestIterations = infoStart.iterations;
        end
    end
end
