%% compare_projected_mra.m
%
% Compare four reconstruction methods for projected MRA:
%
%   1. EM
%   2. Direct optimization over projected T-moments using lsqnonlin
%   3. Direct optimization over Fourier-cosine M-moments using lsqnonlin
%   4. Hard-IRLS linear phase synchronization
%
% This version additionally records:
%   - raw runtime
%   - iteration counts
%   - raw time per iteration
%   - normalized runtime using a reference time-per-iteration estimated
%     from the first few sigma values
%
% Model:
%   y_i = Pi(R_{ell_i} theta) + xi_i,
%   xi_i ~ N(0, sigma^2 I_q),
%   p = 2q+1,
%   Pi(v)[j] = v[j] + v[-j].
%
% Requires:
%   Optimization Toolbox for lsqnonlin.

clear; clc; close all;

% Make all helper functions in this split version available.
thisFile = mfilename('fullpath');
if ~isempty(thisFile)
    addpath(genpath(fileparts(thisFile)));
else
    addpath(genpath(pwd));
end

%% ========================================================================
%                           User parameters
% ========================================================================

rng(1);

% Problem size
p = 13;
q = (p - 1)/2;

% Samples and trials
n = 100000;
numTrials = 100;

% Noise values
sigmaList = logspace(log10(0.05),0,20);

% Success threshold for dihedral-orbit error
successErrTol = 0.20;

% Hard-IRLS settings
irls.numIters = 25;
irls.tol = 1e-8;
irls.minIters = 2;
irls.gaugeWeight = 1e-4;
irls.clipCosines = true;

% Direct moment optimization settings, using lsqnonlin
opt.maxIter = 300;
opt.maxFunEvals = 3000;
opt.display = 'off';
opt.algorithm = 'levenberg-marquardt';

% Number of extra random starts for direct T and direct M optimizations.
% Each method always uses hard-IRLS as its first start.
numRandomMomentStarts = 20;

% EM settings
em.numIters = 2000;
em.tol = 1e-8;
em.numRandomStarts = 5;
em.includeIRLSStart = true;

% Iteration/runtime normalization settings.
% The reference time-per-iteration is estimated from the first few sigma
% points, before long-run slowdown becomes significant.
runtimeCalibNumSigmas = 3;

% Debug printing
debug = true;
debugFirstTrialEachSigma = true;

% Plot/save options
makePlots = true;
savePlots = false;
plotFolder = 'projected_mra_four_method_lsqnonlin_results';

if savePlots && ~exist(plotFolder, 'dir')
    mkdir(plotFolder);
end

if exist('lsqnonlin', 'file') ~= 2
    error('This script requires lsqnonlin from MATLAB Optimization Toolbox.');
end

fprintf('\n============================================================\n');
fprintf('Projected MRA four-method comparison using lsqnonlin\n');
fprintf('p=%d, q=%d, n=%d, trials=%d\n', p, q, n, numTrials);
fprintf('Methods: EM, direct T, direct M, hard IRLS\n');
fprintf('Runtime normalization uses first %d sigma-points.\n', runtimeCalibNumSigmas);
fprintf('============================================================\n\n');

%% ========================================================================
%                         Generate generic signal
% ========================================================================

thetaTrue = generate_generic_real_signal(p, 0.05);
thetaTrue = thetaTrue / norm(thetaTrue);

A = cosine_matrix_A(p);
Ainv = inv(A);

thetaHatTrue = fft(thetaTrue) / sqrt(p);
rTrue = abs(thetaHatTrue(2:q+1));

fprintf('Generated signal. min nonzero Fourier magnitude = %.4g\n\n', min(rTrue));

%% ========================================================================
%                             Storage
% ========================================================================

numSigmas = numel(sigmaList);

% Errors
errIRLS = nan(numSigmas, numTrials);
errDirectT = nan(numSigmas, numTrials);
errDirectM = nan(numSigmas, numTrials);
errEM = nan(numSigmas, numTrials);

% Raw runtimes
timePreprocess = nan(numSigmas, numTrials);
timeIRLS = nan(numSigmas, numTrials);
timeDirectT = nan(numSigmas, numTrials);
timeDirectM = nan(numSigmas, numTrials);
timeEM = nan(numSigmas, numTrials);

% Iteration counts
iterIRLS = nan(numSigmas, numTrials);

iterDirectT_total = nan(numSigmas, numTrials);
iterDirectM_total = nan(numSigmas, numTrials);
iterEM_total = nan(numSigmas, numTrials);

iterDirectT_best = nan(numSigmas, numTrials);
iterDirectM_best = nan(numSigmas, numTrials);
iterEM_best = nan(numSigmas, numTrials);

funcDirectT_total = nan(numSigmas, numTrials);
funcDirectM_total = nan(numSigmas, numTrials);
funcDirectT_best = nan(numSigmas, numTrials);
funcDirectM_best = nan(numSigmas, numTrials);

% Raw time per iteration
timePerIterIRLS = nan(numSigmas, numTrials);
timePerIterDirectT = nan(numSigmas, numTrials);
timePerIterDirectM = nan(numSigmas, numTrials);
timePerIterEM = nan(numSigmas, numTrials);

% Normalized runtimes, filled after the Monte Carlo loop
normTimeIRLS = nan(numSigmas, numTrials);
normTimeDirectT = nan(numSigmas, numTrials);
normTimeDirectM = nan(numSigmas, numTrials);
normTimeEM = nan(numSigmas, numTrials);

% Losses/objectives
lossDirectT = nan(numSigmas, numTrials);
lossDirectM = nan(numSigmas, numTrials);
llEM = nan(numSigmas, numTrials);

% Diagnostics for moment noise
relT3err = nan(numSigmas, numTrials);
relM3err = nan(numSigmas, numTrials);

%% ========================================================================
%                           Monte Carlo loop
% ========================================================================

for sIdx = 1:numSigmas

    sigmaNoise = sigmaList(sIdx);

    fprintf('\n------------------------------------------------------------\n');
    fprintf('sigma = %.4f, nominal SNR = %.3f\n', sigmaNoise, 1/sigmaNoise^2);
    fprintf('------------------------------------------------------------\n');

    % True noiseless moments for diagnostics
    [~, ~, T3true] = projected_moments_from_theta(thetaTrue);
    [~, M3true] = cosine_moments_from_theta(thetaTrue);

    for trial = 1:numTrials

        doDebug = debug && debugFirstTrialEachSigma && trial == 1;

        %% ------------------------------------------------------------
        % Simulate data
        % ------------------------------------------------------------

        Y = simulate_projected_mra(thetaTrue, n, sigmaNoise);

        %% ------------------------------------------------------------
        % Moment preprocessing
        % ------------------------------------------------------------

        tic;
        [T1hat, T2hat, T3hat] = debiased_projected_moments(Y, sigmaNoise);
        [M2hat, M3hat] = transfer_to_cosine_moments(T1hat, T2hat, T3hat, Ainv, p);
        timePreprocess(sIdx, trial) = toc;

        relT3err(sIdx, trial) = norm(T3hat(:) - T3true(:)) / max(norm(T3true(:)), 1e-12);
        relM3err(sIdx, trial) = norm(M3hat(:) - M3true(:)) / max(norm(M3true(:)), 1e-12);

        %% ------------------------------------------------------------
        % 1. Hard-IRLS phase synchronization
        % ------------------------------------------------------------

        tic;
        [thetaIRLS, infoIRLS] = hard_irls_reconstruction(T1hat, M2hat, M3hat, p, irls);
        timeIRLS(sIdx, trial) = toc;

        iterIRLS(sIdx, trial) = infoIRLS.iterations;
        timePerIterIRLS(sIdx, trial) = timeIRLS(sIdx, trial) / max(infoIRLS.iterations, 1);

        errIRLS(sIdx, trial) = dihedral_orbit_error(thetaIRLS, thetaTrue);

        %% ------------------------------------------------------------
        % 2. Direct T-moment optimization using lsqnonlin
        % ------------------------------------------------------------

        tic;
        [thetaT, fvalT, exitflagT, outputT, infoT] = direct_T_moment_optimization_lsqnonlin( ...
            T1hat, T2hat, T3hat, thetaIRLS, numRandomMomentStarts, opt);
        timeDirectT(sIdx, trial) = toc;

        iterDirectT_total(sIdx, trial) = infoT.totalIterations;
        iterDirectT_best(sIdx, trial) = infoT.bestIterations;
        funcDirectT_total(sIdx, trial) = infoT.totalFuncCount;
        funcDirectT_best(sIdx, trial) = infoT.bestFuncCount;
        timePerIterDirectT(sIdx, trial) = timeDirectT(sIdx, trial) / max(infoT.totalIterations, 1);

        errDirectT(sIdx, trial) = dihedral_orbit_error(thetaT, thetaTrue);
        lossDirectT(sIdx, trial) = fvalT;

        %% ------------------------------------------------------------
        % 3. Direct M-moment optimization using lsqnonlin
        % ------------------------------------------------------------

        tic;
        [thetaM, fvalM, exitflagM, outputM, infoM] = direct_M_moment_optimization_lsqnonlin( ...
            T1hat, M2hat, M3hat, thetaIRLS, numRandomMomentStarts, opt);
        timeDirectM(sIdx, trial) = toc;

        iterDirectM_total(sIdx, trial) = infoM.totalIterations;
        iterDirectM_best(sIdx, trial) = infoM.bestIterations;
        funcDirectM_total(sIdx, trial) = infoM.totalFuncCount;
        funcDirectM_best(sIdx, trial) = infoM.bestFuncCount;
        timePerIterDirectM(sIdx, trial) = timeDirectM(sIdx, trial) / max(infoM.totalIterations, 1);

        errDirectM(sIdx, trial) = dihedral_orbit_error(thetaM, thetaTrue);
        lossDirectM(sIdx, trial) = fvalM;

        %% ------------------------------------------------------------
        % 4. EM
        % ------------------------------------------------------------

        tic;
        [thetaEM, bestLL, infoEM] = run_em_with_starts(Y, sigmaNoise, thetaIRLS, em);
        timeEM(sIdx, trial) = toc;

        iterEM_total(sIdx, trial) = infoEM.totalIterations;
        iterEM_best(sIdx, trial) = infoEM.bestIterations;
        timePerIterEM(sIdx, trial) = timeEM(sIdx, trial) / max(infoEM.totalIterations, 1);

        errEM(sIdx, trial) = dihedral_orbit_error(thetaEM, thetaTrue);
        llEM(sIdx, trial) = bestLL;

        %% ------------------------------------------------------------
        % Debug output
        % ------------------------------------------------------------

        if doDebug
            fprintf('\n[debug] sigma=%.3f, trial=%d\n', sigmaNoise, trial);
            fprintf('  preprocessing time: %.4f sec\n', timePreprocess(sIdx, trial));
            fprintf('  rel T3 error: %.4g\n', relT3err(sIdx, trial));
            fprintf('  rel M3 error: %.4g\n', relM3err(sIdx, trial));

            fprintf('  IRLS:     err %.4g, time %.4f sec, iters %d, time/iter %.4g sec, converged %d\n', ...
                errIRLS(sIdx, trial), timeIRLS(sIdx, trial), ...
                iterIRLS(sIdx, trial), timePerIterIRLS(sIdx, trial), ...
                infoIRLS.converged);

            fprintf('  Direct T: err %.4g, loss %.4g, time %.4f sec, exitflag %d, best iters %d, total iters %d, time/iter %.4g sec\n', ...
                errDirectT(sIdx, trial), lossDirectT(sIdx, trial), ...
                timeDirectT(sIdx, trial), exitflagT, ...
                iterDirectT_best(sIdx, trial), iterDirectT_total(sIdx, trial), ...
                timePerIterDirectT(sIdx, trial));

            fprintf('  Direct M: err %.4g, loss %.4g, time %.4f sec, exitflag %d, best iters %d, total iters %d, time/iter %.4g sec\n', ...
                errDirectM(sIdx, trial), lossDirectM(sIdx, trial), ...
                timeDirectM(sIdx, trial), exitflagM, ...
                iterDirectM_best(sIdx, trial), iterDirectM_total(sIdx, trial), ...
                timePerIterDirectM(sIdx, trial));

            fprintf('  EM:       err %.4g, LL %.4g, time %.4f sec, best iters %d, total iters %d, time/iter %.4g sec\n', ...
                errEM(sIdx, trial), llEM(sIdx, trial), timeEM(sIdx, trial), ...
                iterEM_best(sIdx, trial), iterEM_total(sIdx, trial), ...
                timePerIterEM(sIdx, trial));
        end
    end

    %% Per-sigma summary
    fprintf('\nSummary for sigma = %.3f:\n', sigmaNoise);
    fprintf('  median errors: IRLS %.4g | T %.4g | M %.4g | EM %.4g\n', ...
        nanmedian_local(errIRLS(sIdx,:)), ...
        nanmedian_local(errDirectT(sIdx,:)), ...
        nanmedian_local(errDirectM(sIdx,:)), ...
        nanmedian_local(errEM(sIdx,:)));

    fprintf('  median raw times:  IRLS %.4g | T %.4g | M %.4g | EM %.4g | preprocess %.4g sec\n', ...
        nanmedian_local(timeIRLS(sIdx,:)), ...
        nanmedian_local(timeDirectT(sIdx,:)), ...
        nanmedian_local(timeDirectM(sIdx,:)), ...
        nanmedian_local(timeEM(sIdx,:)), ...
        nanmedian_local(timePreprocess(sIdx,:)));

    fprintf('  avg total iters:  IRLS %.3g | T %.3g | M %.3g | EM %.3g\n', ...
        nanmean_local(iterIRLS(sIdx,:)), ...
        nanmean_local(iterDirectT_total(sIdx,:)), ...
        nanmean_local(iterDirectM_total(sIdx,:)), ...
        nanmean_local(iterEM_total(sIdx,:)));
end

%% ========================================================================
%             Runtime normalization by reference time per iteration
% ========================================================================

calibIdx = 1:min(runtimeCalibNumSigmas, numSigmas);

refTimePerIterIRLS = nanmedian_local(timePerIterIRLS(calibIdx,:));
refTimePerIterDirectT = nanmedian_local(timePerIterDirectT(calibIdx,:));
refTimePerIterDirectM = nanmedian_local(timePerIterDirectM(calibIdx,:));
refTimePerIterEM = nanmedian_local(timePerIterEM(calibIdx,:));

normTimeIRLS = iterIRLS * refTimePerIterIRLS;
normTimeDirectT = iterDirectT_total * refTimePerIterDirectT;
normTimeDirectM = iterDirectM_total * refTimePerIterDirectM;
normTimeEM = iterEM_total * refTimePerIterEM;

fprintf('\n============================================================\n');
fprintf('Reference time per iteration, estimated from first %d sigma-points\n', numel(calibIdx));
fprintf('============================================================\n');
fprintf('IRLS:     %.4g sec/iter\n', refTimePerIterIRLS);
fprintf('Direct T: %.4g sec/iter\n', refTimePerIterDirectT);
fprintf('Direct M: %.4g sec/iter\n', refTimePerIterDirectM);
fprintf('EM:       %.4g sec/iter\n', refTimePerIterEM);

%% ========================================================================
%                            Final summaries
% ========================================================================

successIRLS = mean_success(errIRLS, successErrTol);
successT = mean_success(errDirectT, successErrTol);
successM = mean_success(errDirectM, successErrTol);
successEM = mean_success(errEM, successErrTol);

fprintf('\n\n============================================================\n');
fprintf('Final summary: median error\n');
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s\n', 'sigma', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.4g %12.4g %12.4g %12.4g\n', ...
        sigmaList(sIdx), ...
        nanmedian_local(errIRLS(sIdx,:)), ...
        nanmedian_local(errDirectT(sIdx,:)), ...
        nanmedian_local(errDirectM(sIdx,:)), ...
        nanmedian_local(errEM(sIdx,:)));
end

fprintf('\n============================================================\n');
fprintf('Final summary: MSE\n');
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s\n', 'sigma', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.4g %12.4g %12.4g %12.4g\n', ...
        sigmaList(sIdx), ...
        nanmean_local(errIRLS(sIdx,:).^2), ...
        nanmean_local(errDirectT(sIdx,:).^2), ...
        nanmean_local(errDirectM(sIdx,:).^2), ...
        nanmean_local(errEM(sIdx,:).^2));
end

fprintf('\n============================================================\n');
fprintf('Final summary: success probability, error < %.3f\n', successErrTol);
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s\n', 'sigma', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.3f %12.3f %12.3f %12.3f\n', ...
        sigmaList(sIdx), ...
        successIRLS(sIdx), successT(sIdx), successM(sIdx), successEM(sIdx));
end

fprintf('\n============================================================\n');
fprintf('Final summary: median raw running time, seconds\n');
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s %12s\n', ...
    'sigma', 'preproc', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.4g %12.4g %12.4g %12.4g %12.4g\n', ...
        sigmaList(sIdx), ...
        nanmedian_local(timePreprocess(sIdx,:)), ...
        nanmedian_local(timeIRLS(sIdx,:)), ...
        nanmedian_local(timeDirectT(sIdx,:)), ...
        nanmedian_local(timeDirectM(sIdx,:)), ...
        nanmedian_local(timeEM(sIdx,:)));
end

fprintf('\n============================================================\n');
fprintf('Final summary: average total iterations\n');
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s\n', 'sigma', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.3f %12.3f %12.3f %12.3f\n', ...
        sigmaList(sIdx), ...
        nanmean_local(iterIRLS(sIdx,:)), ...
        nanmean_local(iterDirectT_total(sIdx,:)), ...
        nanmean_local(iterDirectM_total(sIdx,:)), ...
        nanmean_local(iterEM_total(sIdx,:)));
end

fprintf('\n============================================================\n');
fprintf('Final summary: median normalized runtime, seconds\n');
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s\n', 'sigma', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.4g %12.4g %12.4g %12.4g\n', ...
        sigmaList(sIdx), ...
        nanmedian_local(normTimeIRLS(sIdx,:)), ...
        nanmedian_local(normTimeDirectT(sIdx,:)), ...
        nanmedian_local(normTimeDirectM(sIdx,:)), ...
        nanmedian_local(normTimeEM(sIdx,:)));
end

fprintf('\n============================================================\n');
fprintf('Final summary: median raw time per iteration, seconds\n');
fprintf('============================================================\n');
fprintf('%8s %12s %12s %12s %12s\n', 'sigma', 'IRLS', 'DirectT', 'DirectM', 'EM');

for sIdx = 1:numSigmas
    fprintf('%8.3f %12.4g %12.4g %12.4g %12.4g\n', ...
        sigmaList(sIdx), ...
        nanmedian_local(timePerIterIRLS(sIdx,:)), ...
        nanmedian_local(timePerIterDirectT(sIdx,:)), ...
        nanmedian_local(timePerIterDirectM(sIdx,:)), ...
        nanmedian_local(timePerIterEM(sIdx,:)));
end

%% ========================================================================
%                                Plots
% ========================================================================

if makePlots

    %% Median error
    figure('Name','Median reconstruction error');
    semilogy(sigmaList, col_nanmedian(errIRLS), '-o', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmedian(errDirectT), '-s', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(errDirectM), '-d', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(errEM), '-p', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel('Median dihedral-orbit error');
    title(sprintf('Reconstruction error, p=%d, q=%d, n=%d', p, q, n));
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'median_error.png'));
    end

    %% MSE
    figure('Name','MSE');
    semilogy(sigmaList, col_nanmean(errIRLS.^2), '-o', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmean(errDirectT.^2), '-s', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmean(errDirectM.^2), '-d', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmean(errEM.^2), '-p', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel('Mean squared dihedral-orbit error');
    title(sprintf('MSE, p=%d, q=%d, n=%d', p, q, n));
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'mse.png'));
    end

    %% Success probability
    figure('Name','Success probability');
    plot(sigmaList, successIRLS, '-o', 'LineWidth', 1.8); hold on;
    plot(sigmaList, successT, '-s', 'LineWidth', 1.8);
    plot(sigmaList, successM, '-d', 'LineWidth', 1.8);
    plot(sigmaList, successEM, '-p', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel(sprintf('Success probability, error < %.2f', successErrTol));
    title(sprintf('Success probability, p=%d, q=%d, n=%d', p, q, n));
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'southwest');
    ylim([0 1.05]);
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'success_probability.png'));
    end

    %% Raw running time
    figure('Name','Median raw runtime');
    semilogy(sigmaList, col_nanmedian(timeIRLS), '-o', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmedian(timeDirectT), '-s', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(timeDirectM), '-d', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(timeEM), '-p', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(timePreprocess), '--k', 'LineWidth', 1.4);
    xlabel('\sigma');
    ylabel('Median raw runtime, seconds');
    title(sprintf('Raw runtime comparison, p=%d, q=%d, n=%d', p, q, n));
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM', 'Moment preprocessing'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'raw_runtime.png'));
    end

    %% Normalized runtime
    figure('Name','Normalized median runtime');
    semilogy(sigmaList, col_nanmedian(normTimeIRLS), '-o', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmedian(normTimeDirectT), '-s', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(normTimeDirectM), '-d', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(normTimeEM), '-p', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel('Median normalized runtime, seconds');
    title(sprintf('Normalized runtime using first %d sigma-points', numel(calibIdx)));
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'normalized_runtime.png'));
    end

    %% Iteration count
    figure('Name','Average iteration count');
    semilogy(sigmaList, col_nanmean(iterIRLS), '-o', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmean(iterDirectT_total), '-s', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmean(iterDirectM_total), '-d', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmean(iterEM_total), '-p', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel('Average total iterations');
    title(sprintf('Effective iteration count, p=%d, q=%d, n=%d', p, q, n));
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'iteration_count.png'));
    end

    %% Raw time per iteration
    figure('Name','Raw time per iteration');
    semilogy(sigmaList, col_nanmedian(timePerIterIRLS), '-o', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmedian(timePerIterDirectT), '-s', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(timePerIterDirectM), '-d', 'LineWidth', 1.8);
    semilogy(sigmaList, col_nanmedian(timePerIterEM), '-p', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel('Median raw time per iteration, seconds');
    title('Measured time per iteration');
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'raw_time_per_iteration.png'));
    end

    %% Accuracy versus normalized runtime
    figure('Name','Error-normalized-time scatter');
    loglog(normTimeIRLS(:), errIRLS(:), 'o'); hold on;
    loglog(normTimeDirectT(:), errDirectT(:), 's');
    loglog(normTimeDirectM(:), errDirectM(:), 'd');
    loglog(normTimeEM(:), errEM(:), 'p');
    xlabel('Normalized runtime, seconds');
    ylabel('Dihedral-orbit error');
    title('Accuracy versus normalized runtime');
    legend({'Hard IRLS', 'Direct T moments', 'Direct M moments', 'EM'}, ...
        'Location', 'southwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'error_normalized_time_scatter.png'));
    end

    %% T3 vs M3 diagnostic
    figure('Name','T3 versus M3 moment error');
    semilogy(sigmaList, col_nanmedian(relT3err), '-s', 'LineWidth', 1.8); hold on;
    semilogy(sigmaList, col_nanmedian(relM3err), '-d', 'LineWidth', 1.8);
    xlabel('\sigma');
    ylabel('Median relative third-moment error');
    title('Noise amplification in Fourier-cosine moment transform');
    legend({'Projected T^{(3)}', 'Fourier-cosine M^{(3)}'}, ...
        'Location', 'northwest');
    grid on;

    if savePlots
        saveas(gcf, fullfile(plotFolder, 'T3_M3_error_diagnostic.png'));
    end
end
