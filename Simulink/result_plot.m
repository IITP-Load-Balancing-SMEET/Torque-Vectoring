% Load data
Fy_sinus_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyfl_sinus_est.mat').ans.time;
Fyfl_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyfl_sinus_est.mat').ans.Data;
Fyfr_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyfr_sinus_est.mat').ans.Data;
Fyrl_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyrl_sinus_est.mat').ans.Data;
Fyrr_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyrr_sinus_est.mat').ans.Data;

My_sinus_est_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\My_sinus_est.mat').ans.time;
My_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\My_sinus_est.mat').ans.Data;

true_sinus_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyfl_sinus_true.mat').ans.time;
Fyfl_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyfl_sinus_true.mat').ans.Data;
Fyfr_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyfr_sinus_true.mat').ans.Data;
Fyrl_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyrl_sinus_true.mat').ans.Data;
Fyrr_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\Fyrr_sinus_true.mat').ans.Data;

My_sinus_true_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\My_sinus_true.mat').ans.time;
My_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\sinus\My_sinus_true.mat').ans.Data;

Fyfl_sinus_true_interp = interp1(true_sinus_time, Fyfl_sinus_true, Fy_sinus_time, 'linear', 'extrap');
Fyfr_sinus_true_interp = interp1(true_sinus_time, Fyfr_sinus_true, Fy_sinus_time, 'linear', 'extrap');
Fyrl_sinus_true_interp = interp1(true_sinus_time, Fyrl_sinus_true, Fy_sinus_time, 'linear', 'extrap');
Fyrr_sinus_true_interp = interp1(true_sinus_time, Fyrr_sinus_true, Fy_sinus_time, 'linear', 'extrap');
My_sinus_true_interp = interp1(My_sinus_true_time, My_sinus_true, My_sinus_est_time, 'linear', 'extrap');

% Compute RMSE
rmse_sinus_fl = sqrt(mean((Fyfl_sinus_est - Fyfl_sinus_true_interp).^2));
rmse_sinus_fr = sqrt(mean((Fyfr_sinus_est - Fyfr_sinus_true_interp).^2));
rmse_sinus_rl = sqrt(mean((Fyrl_sinus_est - Fyrl_sinus_true_interp).^2));
rmse_sinus_rr = sqrt(mean((Fyrr_sinus_est - Fyrr_sinus_true_interp).^2));
rmse_sinus_my = sqrt(mean((My_sinus_est - My_sinus_true_interp).^2));

fontSize = 18;
legendFontSize = 12;
time_limit = 60;
sinus_indices = Fy_sinus_time <= time_limit;
My_sinus_indices = My_sinus_est_time <= time_limit;

% Print the RMSE values
fprintf('Fyfl error: %.2f, Fyfr error: %.2f, Fyrl error: %.2f, Fyrr error: %.2f, My error: %.2f\n', rmse_sinus_fl, rmse_sinus_fr, rmse_sinus_rl, rmse_sinus_rr, rmse_sinus_my);

% Plotting with small margins between subplots
figure(1);
set(gcf, 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Fy_sinus_time(sinus_indices), Fyfl_sinus_true_interp(sinus_indices)); hold on;
plot(Fy_sinus_time(sinus_indices), Fyfl_sinus_est(sinus_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{fl} (N)', 'FontSize', fontSize);
legend('Fy_{fl, Ref}', 'Fy_{fl, Est}', 'FontSize', legendFontSize);

nexttile;
plot(Fy_sinus_time(sinus_indices), Fyfr_sinus_true_interp(sinus_indices)); hold on;
plot(Fy_sinus_time(sinus_indices), Fyfr_sinus_est(sinus_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{fr} (N)', 'FontSize', fontSize);
legend('Fy_{fr, Ref}', 'Fy_{fr, Est}', 'FontSize', legendFontSize);

saveas(figure(1), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\frontForce_sinus.png');

figure(2);
set(gcf, 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Fy_sinus_time(sinus_indices), Fyrl_sinus_true_interp(sinus_indices)); hold on;
plot(Fy_sinus_time(sinus_indices), Fyrl_sinus_est(sinus_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{rl} (N)', 'FontSize', fontSize);
legend('Fy_{rl, Ref}', 'Fy_{rl, Est}', 'FontSize', legendFontSize);

nexttile;
plot(Fy_sinus_time(sinus_indices), Fyrr_sinus_true_interp(sinus_indices)); hold on;
plot(Fy_sinus_time(sinus_indices), Fyrr_sinus_est(sinus_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{rr} (N)', 'FontSize', fontSize);
legend( 'Fy_{rr, Ref}', 'Fy_{rr, Est}', 'FontSize', legendFontSize);

saveas(figure(2), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\rearForce_sinus.png');

figure(3);
set(gcf, 'Position', [100, 100, 1000, 400]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 

nexttile;
plot(My_sinus_est_time(My_sinus_indices), My_sinus_true_interp(My_sinus_indices)); hold on;
plot(My_sinus_est_time(My_sinus_indices), My_sinus_est(My_sinus_indices)); hold off;
xlabel('Time (s)','FontSize', fontSize);
ylabel('M_{y} (N.m)', 'FontSize', fontSize);
legend('M_{y, Ref}', 'M_{y, Est}', 'FontSize', legendFontSize);

saveas(figure(3), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\My_sinus.png');


% Load data
Fy_steady_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyfl_steady_est.mat').ans.time;
Fyfl_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyfl_steady_est.mat').ans.Data;
Fyfr_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyfr_steady_est.mat').ans.Data;
Fyrl_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyrl_steady_est.mat').ans.Data;
Fyrr_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyrr_steady_est.mat').ans.Data;

My_steady_est_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\My_steady_est.mat').ans.time;
My_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\My_steady_est.mat').ans.Data;

true_steady_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyfl_steady_true.mat').ans.time;
Fyfl_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyfl_steady_true.mat').ans.Data;
Fyfr_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyfr_steady_true.mat').ans.Data;
Fyrl_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyrl_steady_true.mat').ans.Data;
Fyrr_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\Fyrr_steady_true.mat').ans.Data;

My_steady_true_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\My_steady_true.mat').ans.time;
My_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\steady\My_steady_true.mat').ans.Data;

% Interpolate the true values to match the time points of the estimated values
Fyfl_steady_true_interp = interp1(true_steady_time, Fyfl_steady_true, Fy_steady_time, 'linear', 'extrap');
Fyfr_steady_true_interp = interp1(true_steady_time, Fyfr_steady_true, Fy_steady_time, 'linear', 'extrap');
Fyrl_steady_true_interp = interp1(true_steady_time, Fyrl_steady_true, Fy_steady_time, 'linear', 'extrap');
Fyrr_steady_true_interp = interp1(true_steady_time, Fyrr_steady_true, Fy_steady_time, 'linear', 'extrap');
My_steady_true_interp = interp1(My_steady_true_time, My_steady_true, My_steady_est_time, 'linear', 'extrap');

% Compute RMSE
rmse_steady_fl = sqrt(mean((Fyfl_steady_est - Fyfl_steady_true_interp).^2));
rmse_steady_fr = sqrt(mean((Fyfr_steady_est - Fyfr_steady_true_interp).^2));
rmse_steady_rl = sqrt(mean((Fyrl_steady_est - Fyrl_steady_true_interp).^2));
rmse_steady_rr = sqrt(mean((Fyrr_steady_est - Fyrr_steady_true_interp).^2));
rmse_steady_my = sqrt(mean((My_steady_est - My_steady_true_interp).^2));

% Print the RMSE values
fprintf('Fyfl error: %.2f, Fyfr error: %.2f, Fyrl error: %.2f, Fyrr error: %.2f, My error: %.2f\n', rmse_steady_fl, rmse_steady_fr, rmse_steady_rl, rmse_steady_rr, rmse_steady_my);

time_limit = 60;
steady_indices = Fy_steady_time <= time_limit;
My_steady_indices = My_steady_est_time <= time_limit;

% Plotting with small margins between subplots
figure(4);
set(gcf, 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Fy_steady_time(steady_indices), Fyfl_steady_true_interp(steady_indices)); hold on;
plot(Fy_steady_time(steady_indices), Fyfl_steady_est(steady_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{fl} (N)', 'FontSize', fontSize);
legend('Fy_{fl, Est}', 'Fy_{fl, Ref}', 'FontSize', legendFontSize);

nexttile;
plot(Fy_steady_time(steady_indices), Fyfr_steady_true_interp(steady_indices)); hold on;
plot(Fy_steady_time(steady_indices), Fyfr_steady_est(steady_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{fr} (N)', 'FontSize', fontSize);
legend('Fy_{fr, Ref}', 'Fy_{fr, Est}', 'FontSize', legendFontSize);

saveas(figure(4), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\frontForce_steady.png');

figure(5);
set(gcf, 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Fy_steady_time(steady_indices), Fyrl_steady_true_interp(steady_indices)); hold on;
plot(Fy_steady_time(steady_indices), Fyrl_steady_est(steady_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{rl} (N)', 'FontSize', fontSize);
legend('Fy_{rl, Ref}', 'Fy_{rl, Est}', 'FontSize', legendFontSize);

nexttile;
plot(Fy_steady_time(steady_indices), Fyrr_steady_true_interp(steady_indices)); hold on;
plot(Fy_steady_time(steady_indices), Fyrr_steady_est(steady_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{rr} (N)', 'FontSize', fontSize);
legend('Fy_{rr, Ref}', 'Fy_{rr, Est}', 'FontSize', legendFontSize);

saveas(figure(5), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\rearForce_steady.png');

figure(6);
set(gcf, 'Position', [100, 100, 1000, 400]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 

nexttile;
plot(My_steady_est_time(My_steady_indices), My_steady_true_interp(My_steady_indices)); hold on;
plot(My_steady_est_time(My_steady_indices), My_steady_est(My_steady_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('M_{y} (N.m)', 'FontSize', fontSize);
legend('M_{y, Ref}', 'M_{y, Est}', 'FontSize', legendFontSize);
ylim([-500, 500]);
saveas(figure(6), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\My_steady.png');


% Load data
Fy_race_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyfl_race_est.mat').ans.time;
Fyfl_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyfl_race_est.mat').ans.Data;
Fyfr_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyfr_race_est.mat').ans.Data;
Fyrl_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyrl_race_est.mat').ans.Data;
Fyrr_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyrr_race_est.mat').ans.Data;

My_race_est_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\My_race_est.mat').ans.time;
My_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\My_race_est.mat').ans.Data;

true_race_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyfl_race_true.mat').ans.time;
Fyfl_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyfl_race_true.mat').ans.Data;
Fyfr_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyfr_race_true.mat').ans.Data;
Fyrl_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyrl_race_true.mat').ans.Data;
Fyrr_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\Fyrr_race_true.mat').ans.Data;

My_race_true_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\My_race_true.mat').ans.time;
My_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\race\My_race_true.mat').ans.Data;

Fyfl_race_true_interp = interp1(true_race_time, Fyfl_race_true, Fy_race_time, 'linear', 'extrap');
Fyfr_race_true_interp = interp1(true_race_time, Fyfr_race_true, Fy_race_time, 'linear', 'extrap');
Fyrl_race_true_interp = interp1(true_race_time, Fyrl_race_true, Fy_race_time, 'linear', 'extrap');
Fyrr_race_true_interp = interp1(true_race_time, Fyrr_race_true, Fy_race_time, 'linear', 'extrap');
My_race_true_interp = interp1(My_race_true_time, My_race_true, My_race_est_time, 'linear', 'extrap');

% Compute RMSE
rmse_race_fl = sqrt(mean((Fyfl_race_est - Fyfl_race_true_interp).^2));
rmse_race_fr = sqrt(mean((Fyfr_race_est - Fyfr_race_true_interp).^2));
rmse_race_rl = sqrt(mean((Fyrl_race_est - Fyrl_race_true_interp).^2));
rmse_race_rr = sqrt(mean((Fyrr_race_est - Fyrr_race_true_interp).^2));
rmse_race_my = sqrt(mean((My_race_est - My_race_true_interp).^2));

time_limit = 60;
race_indices = Fy_race_time <= time_limit;
My_race_indices = My_race_est_time <= time_limit;

% Print the RMSE values
fprintf('Fyfl error: %.2f, Fyfr error: %.2f, Fyrl error: %.2f, Fyrr error: %.2f, My error: %.2f\n', rmse_race_fl, rmse_race_fr, rmse_race_rl, rmse_race_rr, rmse_race_my);

% Plotting with small margins between subplots
figure(7);
set(gcf, 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Fy_race_time(race_indices), Fyfl_race_true_interp(race_indices)); hold on;
plot(Fy_race_time(race_indices), Fyfl_race_est(race_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{fl} (N)', 'FontSize', fontSize);
legend('Fy_{fl, Ref}', 'Fy_{fl, Est}', 'FontSize', legendFontSize);

nexttile;
plot(Fy_race_time(race_indices), Fyfr_race_true_interp(race_indices)); hold on;
plot(Fy_race_time(race_indices), Fyfr_race_est(race_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{fr} (N)', 'FontSize', fontSize);
legend('Fy_{fr, Ref}', 'Fy_{fr, Est}', 'FontSize', legendFontSize);

saveas(figure(7), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\frontForce_race.png');

figure(8);
set(gcf, 'Position', [100, 100, 1000, 800]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

nexttile;
plot(Fy_race_time(race_indices), Fyrl_race_true_interp(race_indices)); hold on;
plot(Fy_race_time(race_indices), Fyrl_race_est(race_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{rl} (N)', 'FontSize', fontSize);
legend('Fy_{rl, Ref}', 'Fy_{rl, Est}', 'FontSize', legendFontSize);

nexttile;
plot(Fy_race_time(race_indices), Fyrr_race_true_interp(race_indices)); hold on;
plot(Fy_race_time(race_indices), Fyrr_race_est(race_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('Fy_{rr} (N)', 'FontSize', fontSize);
legend('Fy_{rr, Ref}', 'Fy_{rr, Est}', 'FontSize', legendFontSize);

saveas(figure(8), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\rearForce_race.png');

figure(9);
set(gcf, 'Position', [100, 100, 1000, 400]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 

nexttile;
plot(My_race_est_time(My_race_indices), My_race_true_interp(My_race_indices)); hold on;
plot(My_race_est_time(My_race_indices), My_race_est(My_race_indices)); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('M_{y} (N.m)', 'FontSize', fontSize);
legend('M_{y, Ref}', 'M_{y, Est}', 'FontSize', legendFontSize);

saveas(figure(9), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\My_race.png');