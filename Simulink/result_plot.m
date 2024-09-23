% Load data
Fy_sinus_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyfl_sinus_est.mat').ans.time;
Fyfl_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyfl_sinus_est.mat').ans.Data;
Fyfr_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyfr_sinus_est.mat').ans.Data;
Fyrl_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyrl_sinus_est.mat').ans.Data;
Fyrr_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyrr_sinus_est.mat').ans.Data;

My_sinus_est_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\My_sinus_est.mat').ans.time;
My_sinus_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\My_sinus_est.mat').ans.Data;

true_sinus_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyfl_sinus_true.mat').ans.time;
Fyfl_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyfl_sinus_true.mat').ans.Data;
Fyfr_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyfr_sinus_true.mat').ans.Data;
Fyrl_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyrl_sinus_true.mat').ans.Data;
Fyrr_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\Fyrr_sinus_true.mat').ans.Data;

My_sinus_true_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\My_sinus_true.mat').ans.time;
My_sinus_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\sinus\My_sinus_true.mat').ans.Data;

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

fontSize = 20;
legendFontSize = 12;
time_limit = 30;
sinus_indices = Fy_sinus_time <= time_limit;
My_sinus_indices = My_sinus_est_time <= time_limit;

% Print the RMSE values
fprintf('Fyfl error: %.2f, Fyfr error: %.2f, Fyrl error: %.2f, Fyrr error: %.2f, My error: %.2f\n', rmse_sinus_fl, rmse_sinus_fr, rmse_sinus_rl, rmse_sinus_rr, rmse_sinus_my);

% Plotting with small margins between subplots
figure(1);
set(gcf, 'Position', [100, 100, 1000, 1000]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax1 = nexttile;
plot(Fy_sinus_time(sinus_indices), Fyfl_sinus_true_interp(sinus_indices), 'LineWidth', 1); hold on; 
plot(Fy_sinus_time(sinus_indices), Fyfl_sinus_est(sinus_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, FL} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax1, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

ax2 = nexttile;
plot(Fy_sinus_time(sinus_indices), Fyfr_sinus_true_interp(sinus_indices), 'LineWidth', 1); hold on;
plot(Fy_sinus_time(sinus_indices), Fyfr_sinus_est(sinus_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y,FR} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax2, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(1), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\frontForce_sinus.png');

figure(2);
set(gcf, 'Position', [100, 100, 1000, 1000]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax3 = nexttile;
plot(Fy_sinus_time(sinus_indices), Fyrl_sinus_true_interp(sinus_indices), 'LineWidth', 1); hold on;
plot(Fy_sinus_time(sinus_indices), Fyrl_sinus_est(sinus_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y,RL} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax3, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

ax4 = nexttile;
plot(Fy_sinus_time(sinus_indices), Fyrr_sinus_true_interp(sinus_indices), 'LineWidth', 1); hold on;
plot(Fy_sinus_time(sinus_indices), Fyrr_sinus_est(sinus_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, RR} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax4, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(2), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\rearForce_sinus.png');

figure(3);
set(gcf, 'Position', [100, 100, 1000, 500]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 

ax5 = nexttile;
plot(My_sinus_est_time(My_sinus_indices), My_sinus_true_interp(My_sinus_indices), 'LineWidth', 1); hold on;
plot(My_sinus_est_time(My_sinus_indices), My_sinus_est(My_sinus_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)','FontSize', fontSize);
ylabel('M_{y} (N.m)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax5, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(3), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\My_sinus.png');


% Load data
Fy_steady_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyfl_steady_est.mat').ans.time;
Fyfl_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyfl_steady_est.mat').ans.Data;
Fyfr_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyfr_steady_est.mat').ans.Data;
Fyrl_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyrl_steady_est.mat').ans.Data;
Fyrr_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyrr_steady_est.mat').ans.Data;

My_steady_est_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\My_steady_est.mat').ans.time;
My_steady_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\My_steady_est.mat').ans.Data;

true_steady_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyfl_steady_true.mat').ans.time;
Fyfl_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyfl_steady_true.mat').ans.Data;
Fyfr_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyfr_steady_true.mat').ans.Data;
Fyrl_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyrl_steady_true.mat').ans.Data;
Fyrr_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\Fyrr_steady_true.mat').ans.Data;

My_steady_true_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\My_steady_true.mat').ans.time;
My_steady_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\steady\My_steady_true.mat').ans.Data;

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

steady_indices = Fy_steady_time <= time_limit;
My_steady_indices = My_steady_est_time <= time_limit;

% Plotting with small margins between subplots
figure(4);
set(gcf, 'Position', [100, 100, 1000, 1000]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax6 = nexttile;
plot(Fy_steady_time(steady_indices), Fyfl_steady_true_interp(steady_indices), 'LineWidth', 1); hold on;
plot(Fy_steady_time(steady_indices), Fyfl_steady_est(steady_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, FL} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax6, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

ax7 = nexttile;
plot(Fy_steady_time(steady_indices), Fyfr_steady_true_interp(steady_indices), 'LineWidth', 1); hold on;
plot(Fy_steady_time(steady_indices), Fyfr_steady_est(steady_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, FR} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax7, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(4), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\frontForce_steady.png');

figure(5);
set(gcf, 'Position', [100, 100, 1000, 1000]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax8 = nexttile;
plot(Fy_steady_time(steady_indices), Fyrl_steady_true_interp(steady_indices), 'LineWidth', 1); hold on;
plot(Fy_steady_time(steady_indices), Fyrl_steady_est(steady_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, RL} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax8, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

ax9 = nexttile;
plot(Fy_steady_time(steady_indices), Fyrr_steady_true_interp(steady_indices), 'LineWidth', 1); hold on;
plot(Fy_steady_time(steady_indices), Fyrr_steady_est(steady_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, RR} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax9, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(5), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\rearForce_steady.png');

figure(6);
set(gcf, 'Position', [100, 100, 1000, 500]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 

ax10 = nexttile;
plot(My_steady_est_time(My_steady_indices), My_steady_true_interp(My_steady_indices), 'LineWidth', 1); hold on;
plot(My_steady_est_time(My_steady_indices), My_steady_est(My_steady_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('M_{y} (N.m)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax10, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks
ylim([-400, 400]);
saveas(figure(6), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\My_steady.png');


% Load data
Fy_race_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyfl_race_est.mat').ans.time;
Fyfl_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyfl_race_est.mat').ans.Data;
Fyfr_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyfr_race_est.mat').ans.Data;
Fyrl_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyrl_race_est.mat').ans.Data;
Fyrr_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyrr_race_est.mat').ans.Data;

My_race_est_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\My_race_est.mat').ans.time;
My_race_est = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\My_race_est.mat').ans.Data;

true_race_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyfl_race_true.mat').ans.time;
Fyfl_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyfl_race_true.mat').ans.Data;
Fyfr_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyfr_race_true.mat').ans.Data;
Fyrl_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyrl_race_true.mat').ans.Data;
Fyrr_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\Fyrr_race_true.mat').ans.Data;

My_race_true_time = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\My_race_true.mat').ans.time;
My_race_true = load('C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Simulink\results\AEKF\race\My_race_true.mat').ans.Data;

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
set(gcf, 'Position', [100, 100, 1000, 1000]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax11 = nexttile;
plot(Fy_race_time(race_indices), Fyfl_race_true_interp(race_indices), 'LineWidth', 1); hold on;
plot(Fy_race_time(race_indices), Fyfl_race_est(race_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, FL} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax11, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

ax12 = nexttile;
plot(Fy_race_time(race_indices), Fyfr_race_true_interp(race_indices), 'LineWidth', 1); hold on;
plot(Fy_race_time(race_indices), Fyfr_race_est(race_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, FR} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax12, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(7), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\frontForce_race.png');

figure(8);
set(gcf, 'Position', [100, 100, 1000, 1000]);
tiledlayout(2, 1, 'Padding', 'compact', 'TileSpacing', 'compact');

ax13 = nexttile;
plot(Fy_race_time(race_indices), Fyrl_race_true_interp(race_indices),'LineWidth', 1); hold on;
plot(Fy_race_time(race_indices), Fyrl_race_est(race_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, RL} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax13, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

ax14 = nexttile;
plot(Fy_race_time(race_indices), Fyrr_race_true_interp(race_indices), 'LineWidth', 1); hold on;
plot(Fy_race_time(race_indices), Fyrr_race_est(race_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('F_{y, RR} (N)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax14, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(8), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\rearForce_race.png');

figure(9);
set(gcf, 'Position', [100, 100, 1000, 500]);
tiledlayout(1, 1, 'Padding', 'compact', 'TileSpacing', 'compact'); 

ax15 = nexttile;
plot(My_race_est_time(My_race_indices), My_race_true_interp(My_race_indices), 'LineWidth', 1); hold on;
plot(My_race_est_time(My_race_indices), My_race_est(My_race_indices), 'LineWidth', 1); hold off;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('M_{y} (N.m)', 'FontSize', fontSize);
legend('GroundTruth', 'Estimation', 'FontSize', legendFontSize);
set(ax15, 'FontSize', 15);  % Set font size for all text in the axes, including x-ticks

saveas(figure(9), 'C:\Users\SMEET_SIMUL\Desktop\Torque-Vectoring\Paper\Figures\My_race.png');