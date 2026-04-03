function [BestX, BestF, HisBestFit, out] = WECO(obj_fun, D, N0, MaxFEs, lb, ub, varargin)
% WECO  Water-Energy Complementary Optimizer
%
%   [BestX, BestF, HisBestFit, out] = WECO(obj_fun, D, N0, MaxFEs, lb, ub)
%   minimizes the objective function obj_fun over a D-dimensional bounded
%   search space.
%
%   [BestX, BestF, HisBestFit, out] = WECO(..., opts)
%   accepts an optional options struct with fields:
%       - history_length : number of stored history points
%       - seed           : RNG seed
%       - verbose        : true / false
%
%   [BestX, BestF, HisBestFit, out] = WECO(..., 'seed', 2026, 'verbose', true)
%   also supports traditional name-value pairs.
%
% Inputs
%   obj_fun  : objective function handle. It may accept either a row vector
%              (1-by-D) or a column vector (D-by-1).
%   D        : decision-space dimension.
%   N0       : initial population size.
%   MaxFEs   : maximum number of function evaluations.
%   lb, ub   : lower/upper bounds, each either scalar or 1-by-D vector.
%
% Outputs
%   BestX       : best solution found.
%   BestF       : best objective value found.
%   HisBestFit  : best-so-far history sampled uniformly in evaluation space.
%   out         : struct containing diagnostic information.
%
% Notes
%   1) The algorithm always respects the function-evaluation budget MaxFEs.
%   2) History is logged in evaluation space rather than iteration count.
%   3) The implementation is self-contained and uses only MATLAB built-ins.
%
% Author-prepared open-source version for repository release.

    % Basic validation (kept compatible with older MATLAB releases)
    if nargin < 6
        error('WECO requires at least 6 inputs: obj_fun, D, N0, MaxFEs, lb, ub.');
    end
    if ~isa(obj_fun, 'function_handle')
        error('obj_fun must be a function handle.');
    end
    if ~isscalar(D) || D <= 0 || floor(D) ~= D
        error('D must be a positive integer scalar.');
    end
    if ~isscalar(N0) || N0 <= 3 || floor(N0) ~= N0
        error('N0 must be an integer scalar greater than 3.');
    end
    if ~isscalar(MaxFEs) || MaxFEs <= 0 || floor(MaxFEs) ~= MaxFEs
        error('MaxFEs must be a positive integer scalar.');
    end

    opts = struct();
    opts.history_length = ceil(MaxFEs / N0);
    opts.seed = [];
    opts.verbose = false;

    if ~isempty(varargin)
        if numel(varargin) == 1 && isstruct(varargin{1})
            user_opts = varargin{1};
            fn = fieldnames(user_opts);
            for k = 1:numel(fn)
                opts.(fn{k}) = user_opts.(fn{k});
            end
        else
            if mod(numel(varargin), 2) ~= 0
                error('Optional arguments must be provided as a struct or as name-value pairs.');
            end
            for k = 1:2:numel(varargin)
                name = varargin{k};
                value = varargin{k+1};
                if isstring(name), name = char(name); end
                if ~ischar(name)
                    error('Name-value option names must be character vectors or strings.');
                end
                opts.(name) = value;
            end
        end
    end

    if ~isscalar(opts.history_length) || opts.history_length <= 0
        error('opts.history_length must be a positive integer scalar.');
    end
    opts.history_length = max(1, round(opts.history_length));
    if ~isempty(opts.seed)
        rng(opts.seed);
    end
    opts.verbose = logical(opts.verbose);

    lb = expand_bounds(lb, D, 'lb');
    ub = expand_bounds(ub, D, 'ub');

    if any(ub <= lb)
        error('Each upper bound must be strictly greater than the corresponding lower bound.');
    end

    history_length = max(1, opts.history_length);
    history_eval_points = unique(round(linspace(1, MaxFEs, history_length)));
    history_length = numel(history_eval_points);
    HisBestFit = nan(history_length, 1);
    next_hist_idx = 1;

    % ---------------------------------------------------------------------
    % 0) Initialization and nondimensional scale
    % ---------------------------------------------------------------------
    evals = 0;
    L = norm(ub - lb) / sqrt(D);
    L = max(L, eps);

    X = lb + rand(N0, D) .* (ub - lb);
    F = zeros(N0, 1);
    for i = 1:N0
        F(i) = safe_eval(obj_fun, X(i, :));
    end
    evals = evals + N0;

    [F, idx] = sort(F, 'ascend');
    X = X(idx, :);

    BestF = F(1);
    BestX = X(1, :);
    [HisBestFit, next_hist_idx] = fill_history(HisBestFit, history_eval_points, next_hist_idx, evals, BestF);

    % ---------------------------------------------------------------------
    % 1) Hydrological states and adaptive memories
    % ---------------------------------------------------------------------
    % Soil-moisture state W in [0, 1]
    W = 0.3 + 0.4 * rand(N0, 1);
    W = clip(W, 0, 1);

    % Runoff-response state Q in [0, 1]
    Q = zeros(N0, 1);

    % Success-history memories
    Hmem = 10;
    M_eta = 0.5 * ones(Hmem, 1);   % evolution-intensity memory
    M_kappa = 0.5 * ones(Hmem, 1); % connectivity memory
    mem_idx = 1;

    % Historical retention pool (archive)
    archive.X = zeros(0, D);
    archive.F = zeros(0, 1);
    archive.max_size = N0;

    % Physical baseline parameters
    c_base_raw = 0.01;
    c_base = c_base_raw / L;
    k_c = 2.0;

    NPmin = max(4, round(0.15 * N0));
    sigmaQ_prev = 0.25;

    % Diagnostics
    NP_hist = [];
    Wmean_hist = [];
    Qmean_hist = [];
    sigmaQ_hist = [];
    lambda_hist = [];
    tau_hist = [];

    gen = 0;
    while evals < MaxFEs
        gen = gen + 1;
        nPop = size(X, 1);
        tau = min(1, evals / MaxFEs);

        % Stage-dependent infiltration baseline and heterogeneity index
        c_tau = c_base * (1 - tau) + 1e-5 / L;
        n_tau = clip(0.1 + 9.9 * (tau ^ 2), 0.1, 10.0);
        lambda_gen = min(0.5, (1 - tau) * sigmaQ_prev);

        if opts.verbose
            fprintf('WECO | generation %d | evals %d/%d | best %.6e\n', gen, evals, MaxFEs, BestF);
        end

        X_old = X;
        F_old = F;

        S_eta = [];
        S_kappa = [];
        S_df = [];
        Q_list = nan(nPop, 1);

        for i = 1:nPop
            if evals >= MaxFEs
                break;
            end

            % -------------------------------------------------------------
            % A) Equivalent rainfall input from inter-state difference
            % -------------------------------------------------------------
            r = randperm(nPop, 2);
            while any(r == i)
                r = randperm(nPop, 2);
            end
            R_i = (norm(X_old(r(1), :) - X_old(r(2), :)) / sqrt(D)) / L + eps;

            % -------------------------------------------------------------
            % B) Infiltration capacity and actual infiltration
            % C_i^inf = c(tau) * (1 + k_c * (1 - W_i))
            % A_i^inf = R_i * C_i^inf / (R_i^n + (C_i^inf)^n)^(1/n)
            % -------------------------------------------------------------
            C_i_inf = c_tau * (1 + k_c * (1 - W(i)));
            C_i_inf = max(C_i_inf, 1e-12);
            denom = (R_i ^ n_tau + C_i_inf ^ n_tau) ^ (1 / n_tau);
            A_i_inf = (R_i * C_i_inf) / (denom + eps);

            alpha_i = clip(A_i_inf / (R_i + eps), 0, 1);
            rho_i = clip(1 - alpha_i, 0, 1);

            % -------------------------------------------------------------
            % C) Runoff routing: linear-reservoir style response
            % -------------------------------------------------------------
            Q(i) = clip((1 - tau) * Q(i) + tau * rho_i, 0, 1);
            Q_list(i) = Q(i);

            % -------------------------------------------------------------
            % D) State-dependent control variables
            % -------------------------------------------------------------
            p_i = clip(0.05 + 0.45 * W(i), 0.05, 0.50);
            p_num = max(2, round(p_i * nPop));
            pbest = randi(p_num);

            psi_i = Q(i) - W(i);
            rmem = randi(Hmem);

            mu_eta = clip(M_eta(rmem) + lambda_gen * psi_i, 0.05, 0.95);
            mu_kappa = clip(M_kappa(rmem) - lambda_gen * psi_i, 0.0, 1.0);

            eta_i = cauchy_rand(mu_eta, 0.1);
            retry = 0;
            while eta_i <= 0 && retry < 5
                eta_i = cauchy_rand(mu_eta, 0.1);
                retry = retry + 1;
            end
            if eta_i <= 0
                eta_i = max(0.05, mu_eta);
            end
            eta_i = min(eta_i, 1.0);

            kappa_i = clip(mu_kappa + 0.1 * randn, 0.0, 1.0);

            % -------------------------------------------------------------
            % E) Select tributary and retention-pool states
            % -------------------------------------------------------------
            rr1 = randi(nPop);
            while rr1 == i || rr1 == pbest
                rr1 = randi(nPop);
            end

            use_archive = (~isempty(archive.X)) && (rand < Q(i));
            if use_archive
                rr2 = randi(size(archive.X, 1));
                x_r2 = archive.X(rr2, :);
            else
                rr2 = randi(nPop);
                while rr2 == i || rr2 == pbest || rr2 == rr1
                    rr2 = randi(nPop);
                end
                x_r2 = X_old(rr2, :);
            end
            x_r1 = X_old(rr1, :);
            x_main = X_old(pbest, :);

            % -------------------------------------------------------------
            % F) Storage-support feedback and candidate generation
            % -------------------------------------------------------------
            b_i = max(0, W(i) - Q(i));
            x_tilde = X_old(i, :) + eta_i * (x_main - X_old(i, :)) ...
                                  + eta_i * (x_r1 - x_r2) ...
                                  + b_i * eta_i * (BestX - X_old(i, :));

            jrand = randi(D);
            x_hat = X_old(i, :);
            for d = 1:D
                if rand <= kappa_i || d == jrand
                    x_hat(d) = x_tilde(d);
                end
            end
            x_hat = repair_bounds(x_hat, lb, ub);

            % -------------------------------------------------------------
            % G) Evaluation and greedy acceptance
            % -------------------------------------------------------------
            f_hat = safe_eval(obj_fun, x_hat);
            evals = evals + 1;

            if f_hat < F(i)
                archive = archive_add(archive, X_old(i, :), F(i));
                X(i, :) = x_hat;
                F(i) = f_hat;

                delta_f = abs(F_old(i) - f_hat);
                S_eta(end + 1, 1) = eta_i; %#ok<AGROW>
                S_kappa(end + 1, 1) = kappa_i; %#ok<AGROW>
                S_df(end + 1, 1) = delta_f; %#ok<AGROW>

                if f_hat < BestF
                    BestF = f_hat;
                    BestX = x_hat;
                end
            end

            % -------------------------------------------------------------
            % H) Soil-moisture water-balance closure
            % W_i <- W_i + (1-W_i)*alpha_i - tau*W_i - Q_i*W_i
            % -------------------------------------------------------------
            W(i) = clip(W(i) + (1 - W(i)) * alpha_i - tau * W(i) - Q(i) * W(i), 0, 1);

            [HisBestFit, next_hist_idx] = fill_history(HisBestFit, history_eval_points, next_hist_idx, evals, BestF);
            if evals >= MaxFEs
                break;
            end
        end

        [F, idx] = sort(F, 'ascend');
        X = X(idx, :);
        W = W(idx);
        Q = Q(idx);

        if F(1) < BestF
            BestF = F(1);
            BestX = X(1, :);
        end

        % Update success-history memories
        if ~isempty(S_eta) && sum(S_df) > 0
            weights = S_df ./ (sum(S_df) + eps);
            M_eta(mem_idx) = clip(sum(weights .* (S_eta .^ 2)) / (sum(weights .* S_eta) + eps), 0.05, 0.95);
            M_kappa(mem_idx) = clip(sum(weights .* S_kappa), 0.0, 1.0);
            mem_idx = mem_idx + 1;
            if mem_idx > Hmem
                mem_idx = 1;
            end
        end

        % Update runoff heterogeneity
        Q_valid = Q_list(~isnan(Q_list));
        if numel(Q_valid) >= 2
            sigmaQ_prev = clip(std(Q_valid), 0, 1);
        end

        % Population-size reduction
        tau = min(1, evals / MaxFEs);
        nPop_target = round(NPmin + (N0 - NPmin) * (1 - tau));
        nPop_target = max(NPmin, min(N0, nPop_target));
        if size(X, 1) > nPop_target
            X = X(1:nPop_target, :);
            F = F(1:nPop_target);
            W = W(1:nPop_target);
            Q = Q(1:nPop_target);
        end

        archive.max_size = size(X, 1);
        archive = archive_trim(archive);

        NP_hist(end + 1, 1) = size(X, 1); %#ok<AGROW>
        Wmean_hist(end + 1, 1) = mean(W); %#ok<AGROW>
        Qmean_hist(end + 1, 1) = mean(Q); %#ok<AGROW>
        sigmaQ_hist(end + 1, 1) = sigmaQ_prev; %#ok<AGROW>
        lambda_hist(end + 1, 1) = lambda_gen; %#ok<AGROW>
        tau_hist(end + 1, 1) = tau; %#ok<AGROW>
    end

    if any(isnan(HisBestFit))
        last_valid = find(~isnan(HisBestFit), 1, 'last');
        if isempty(last_valid)
            HisBestFit(:) = BestF;
        else
            HisBestFit(last_valid + 1:end) = BestF;
        end
    end

    out = struct();
    out.EvalsUsed = evals;
    out.MaxFEs = MaxFEs;
    out.HisBestFit = HisBestFit;
    out.Diversity = compute_diversity(X, BestX, lb, ub);

    out.InitialPopulationSize = N0;
    out.FinalPopulationSize = size(X, 1);
    out.ArchiveSize = size(archive.X, 1);

    out.L_scale = L;
    out.c_base_raw = c_base_raw;
    out.c_base_scaled = c_base;
    out.k_c = k_c;
    out.Hmem = Hmem;
    out.M_eta = M_eta;
    out.M_kappa = M_kappa;

    out.W_mean_final = mean(W);
    out.Q_mean_final = mean(Q);
    out.NP_hist = NP_hist;
    out.Wmean_hist = Wmean_hist;
    out.Qmean_hist = Qmean_hist;
    out.sigmaQ_hist = sigmaQ_hist;
    out.lambda_hist = lambda_hist;
    out.tau_hist = tau_hist;
end

% -------------------------------------------------------------------------
% Helpers
% -------------------------------------------------------------------------
function x = expand_bounds(x, D, name)
    if isscalar(x)
        x = repmat(x, 1, D);
    else
        x = x(:).';
    end
    if numel(x) ~= D
        error('Bound %s must be scalar or length-D vector (D = %d).', name, D);
    end
end

function y = safe_eval(obj_fun, x)
    try
        y = obj_fun(x);
    catch
        y = obj_fun(x');
    end

    if ~isscalar(y) || ~isfinite(y)
        if isnan(y) || isinf(y)
            y = 1e30;
        else
            error('Objective function must return one finite scalar value.');
        end
    end

    y = min(max(double(y), -1e30), 1e30);
end

function x = repair_bounds(x, lb, ub)
    above = x > ub;
    below = x < lb;
    if any(above)
        x(above) = lb(above) + rand(1, sum(above)) .* (ub(above) - lb(above));
    end
    if any(below)
        x(below) = lb(below) + rand(1, sum(below)) .* (ub(below) - lb(below));
    end
end

function archive = archive_add(archive, x_old, f_old)
    archive.X(end + 1, :) = x_old;
    archive.F(end + 1, 1) = f_old;
    archive = archive_trim(archive);
end

function archive = archive_trim(archive)
    while size(archive.X, 1) > archive.max_size
        k = randi(size(archive.X, 1));
        archive.X(k, :) = [];
        archive.F(k, :) = [];
    end
end

function y = cauchy_rand(loc, scale)
    y = loc + scale * tan(pi * (rand - 0.5));
end

function y = clip(x, lo, hi)
    y = min(max(x, lo), hi);
end

function diversity = compute_diversity(X, BestX, lb, ub)
    nPop = size(X, 1);
    if nPop <= 1
        diversity = 0;
        return;
    end

    center = mean(X, 1);
    denom = max(norm(ub - lb), eps);
    dist_center = vecnorm(X - center, 2, 2) / denom;
    dist_best = vecnorm(X - BestX, 2, 2) / denom;
    diversity = clip(0.5 * mean(dist_center) + 0.5 * mean(dist_best), 0, 1);
end

function [history, next_idx] = fill_history(history, eval_points, next_idx, evals, bestf)
    while next_idx <= numel(eval_points) && evals >= eval_points(next_idx)
        history(next_idx) = bestf;
        next_idx = next_idx + 1;
    end
end
