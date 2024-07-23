% Load data
Ca0 = 30000;

alpha_fl = load("alpha_fl.mat").ans.Data;
alpha_fr = load("alpha_fr.mat").ans.Data;
alpha_rl = load("alpha_rl.mat").ans.Data;
alpha_rr = load("alpha_rr.mat").ans.Data;

alpha_fl_sinus = load("alpha_fl_sinus.mat").ans.Data;
alpha_fr_sinus = load("alpha_fr_sinus.mat").ans.Data;
alpha_rl_sinus = load("alpha_rl_sinus.mat").ans.Data;
alpha_rr_sinus = load("alpha_rr_sinus.mat").ans.Data;

alpha_fl_sc = load("alpha_fl_sc.mat").ans.Data;
alpha_fr_sc = load("alpha_fr_sc.mat").ans.Data;
alpha_rl_sc = load("alpha_rl_sc.mat").ans.Data;
alpha_rr_sc = load("alpha_rr_sc.mat").ans.Data;

alpha_fl = [alpha_fl; alpha_fl_sinus; alpha_fl_sc];
alpha_fr = [alpha_fr; alpha_fr_sinus; alpha_fr_sc];
alpha_rl = [alpha_rl; alpha_rl_sinus; alpha_rl_sc];
alpha_rr = [alpha_rr; alpha_rr_sinus; alpha_rr_sc];

Fzfl = load("Fzfl.mat").ans.Data;
Fzfr = load("Fzfr.mat").ans.Data;
Fzrl = load("Fzrl.mat").ans.Data;
Fzrr = load("Fzrr.mat").ans.Data;

Fzfl_sinus = load("Fzfl_sinus.mat").ans.Data;
Fzfr_sinus = load("Fzfr_sinus.mat").ans.Data;
Fzrl_sinus = load("Fzrl_sinus.mat").ans.Data;
Fzrr_sinus = load("Fzrr_sinus.mat").ans.Data;

Fzfl_sc = load("Fzfl_sc.mat").ans.Data;
Fzfr_sc = load("Fzfr_sc.mat").ans.Data;
Fzrl_sc = load("Fzrl_sc.mat").ans.Data;
Fzrr_sc = load("Fzrr_sc.mat").ans.Data;

Fzfl = [Fzfl; Fzfl_sinus; Fzfl_sc];
Fzfr = [Fzfr; Fzfr_sinus; Fzfr_sc];
Fzrl = [Fzrl; Fzrl_sinus; Fzrl_sc];
Fzrr = [Fzrr; Fzrr_sinus; Fzrr_sc];

Fyfl = load("Fyfl_2.mat").ans.Data;
Fyfr = load("Fyfr_2.mat").ans.Data;
Fyrl = load("Fyrl_2.mat").ans.Data;
Fyrr = load("Fyrr_2.mat").ans.Data;

Fyfl_sinus = load("Fyfl_2_sinus.mat").ans.Data;
Fyfr_sinus = load("Fyfr_2_sinus.mat").ans.Data;
Fyrl_sinus = load("Fyrl_2_sinus.mat").ans.Data;
Fyrr_sinus = load("Fyrr_2_sinus.mat").ans.Data;

Fyfl_sc = load("Fyfl_2_sc.mat").ans.Data;
Fyfr_sc = load("Fyfr_2_sc.mat").ans.Data;
Fyrl_sc = load("Fyrl_2_sc.mat").ans.Data;
Fyrr_sc = load("Fyrr_2_sc.mat").ans.Data;

Fyfl = [Fyfl; Fyfl_sinus; Fyfl_sc];
Fyfr = [Fyfr; Fyfr_sinus; Fyfr_sc];
Fyrl = [Fyrl; Fyrl_sinus; Fyrl_sc];
Fyrr = [Fyrr; Fyrr_sinus; Fyrr_sc];

data = [alpha_fl, Fzfl, Fyfl; 
        alpha_fr, Fzfr, Fyfr; 
        alpha_rl, Fzrl, Fyrl;
        alpha_rr, Fzrr, Fyrr];

alpha = data(:, 1);
Fz = data(:, 2);
Fy = data(:, 3);


% Set options for lsqnonlin
loss_function = @(params) norm(Fy - -((params(1) * Fz.^2) + (params(2) * Fz) + Ca0) .* alpha);
residual_function = @(params) Fy - -((params(1) * Fz.^2) + (params(2) * Fz) + Ca0) .* alpha;

initial_guess = [1; 1];

bfgs_options = optimoptions('fminunc', 'MaxIterations', 10000, 'MaxFunctionEvaluations', 10000, ...
    "HessianApproximation", 'bfgs', "OptimalityTolerance", 1.0000e-10);
lsq_options = optimoptions('lsqnonlin', 'MaxIterations', 10000, 'MaxFunctionEvaluations', 10000, ...
    'OptimalityTolerance', 1.0000e-10, 'Display', 'iter');

% Perform optimization
[bfgs_optimal_params, fval] = fminunc(loss_function, initial_guess, bfgs_options);
[lsq_optimal_params, resnorm] = lsqnonlin(residual_function, initial_guess, [], [], lsq_options);

% Display the results
bfgs_a_opt = bfgs_optimal_params(1);
bfgs_b_opt = bfgs_optimal_params(2);

lsq_a_opt = lsq_optimal_params(1);
lsq_b_opt = lsq_optimal_params(2);

fprintf('BFGS Optimal parameters: a = %.6f, b = %.6f\n', bfgs_a_opt, bfgs_b_opt);
fprintf('BFGS cost function value: %.4f\n', fval);

fprintf('LSQ Optimal parameters: a = %.6f, b = %.6f\n', lsq_a_opt, lsq_b_opt);
fprintf('LSQ Final residual norm: %.4f\n', resnorm);

% Calculate the estimated forces using the optimized parameters
Fyfl2 = -((optimal_params(1) .* Fzfl.^2) + (optimal_params(2) .* Fzfl) + Ca0) .* alpha_fl;
Fyfr2 = -((optimal_params(1) .* Fzfr.^2) + (optimal_params(2) .* Fzfr) + Ca0) .* alpha_fr;
Fyrl2 = -((optimal_params(1) .* Fzrl.^2) + (optimal_params(2) .* Fzrl) + Ca0) .* alpha_rl;
Fyrr2 = -((optimal_params(1) .* Fzrr.^2) + (optimal_params(2) .* Fzrr) + Ca0) .* alpha_rr;

% Plot the results
figure(1);
subplot(2, 2, 1);
scatter(alpha_fl, Fyfl, 3); hold on;
scatter(alpha_fl, Fyfl2, 3); hold off;
xlabel('Alpha (fl)');
ylabel('Fy (fl)');
title('Front Left');
legend("Fy", "Fy_{est}");

subplot(2, 2, 2);
scatter(alpha_fr, Fyfr, 3); hold on;
scatter(alpha_fr, Fyfr2, 3); hold off;
xlabel('Alpha (fr)');
ylabel('Fy (fr)');
title('Front Right');
legend("Fy", "Fy_{est}");

subplot(2, 2, 3);
scatter(alpha_rl, Fyrl, 3); hold on;
scatter(alpha_rl, Fyrl2, 3); hold off;
xlabel('Alpha (rl)');
ylabel('Fy (rl)');
title('Rear Left');
legend("Fy", "Fy_{est}");

subplot(2, 2, 4);
scatter(alpha_rr, Fyrr, 3); hold on;
scatter(alpha_rr, Fyrr2, 3); hold off;
xlabel('Alpha (rr)');
ylabel('Fy (rr)');
title('Rear Right');
legend("Fy", "Fy_{est}");

figure(2);
subplot(2, 2, 1);
scatter(Fzfl, Fyfl, 3, 'r');
xlabel('Fz (fl)');
ylabel('Fy (fl)');
title('Front Left');

subplot(2, 2, 2);
scatter(Fzfr, Fyfr, 3, 'b');
xlabel('Fz (fr)');
ylabel('Fy (fr)');
title('Front Right');

subplot(2, 2, 3);
scatter(Fzrl, Fyrl, 3, 'g');
xlabel('Fz (rl)');
ylabel('Fy (rl)');
title('Rear Left');

subplot(2, 2, 4);
scatter(Fzrr, Fyrr, 3, 'k');
xlabel('Fz (rr)');
ylabel('Fy (rr)');
title('Rear Right');
