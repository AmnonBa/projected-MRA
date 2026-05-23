function [thetaBest, fBest, exitflagBest, outputBest, info] = direct_T_moment_optimization_lsqnonlin( ...
    T1hat, T2hat, T3hat, thetaInit, numRandomStarts, opt)

    p = numel(thetaInit);

    starts = zeros(p, 1 + numRandomStarts);
    starts(:,1) = thetaInit(:);

    baseNorm = max(norm(thetaInit), 1);
    for s = 1:numRandomStarts
        z = randn(p,1);
        z = z / norm(z) * baseNorm;
        starts(:,1+s) = z;
    end

    options = optimoptions('lsqnonlin', ...
        'Display', opt.display, ...
        'Algorithm', opt.algorithm, ...
        'MaxIterations', opt.maxIter, ...
        'MaxFunctionEvaluations', opt.maxFunEvals, ...
        'FunctionTolerance', 1e-10, ...
        'StepTolerance', 1e-10);

    fBest = inf;
    thetaBest = thetaInit(:);
    exitflagBest = NaN;
    outputBest.iterations = NaN;
    outputBest.funcCount = NaN;

    info.totalIterations = 0;
    info.totalFuncCount = 0;
    info.bestIterations = NaN;
    info.bestFuncCount = NaN;
    info.numSuccessfulStarts = 0;
    info.numStarts = size(starts,2);

    for s = 1:size(starts,2)
        x0 = starts(:,s);
        resfun = @(x) T_moment_residual_vector(x(:), T1hat, T2hat, T3hat);

        try
            [thetaCand, resnorm, ~, exitflag, output] = lsqnonlin(resfun, x0, [], [], options);
        catch ME
            warning('lsqnonlin failed in direct T optimization, start %d: %s', s, ME.message);
            continue;
        end

        info.numSuccessfulStarts = info.numSuccessfulStarts + 1;
        info.totalIterations = info.totalIterations + output.iterations;
        info.totalFuncCount = info.totalFuncCount + output.funcCount;

        if resnorm < fBest
            fBest = resnorm;
            thetaBest = thetaCand(:);
            exitflagBest = exitflag;
            outputBest = output;
            info.bestIterations = output.iterations;
            info.bestFuncCount = output.funcCount;
        end
    end
end
