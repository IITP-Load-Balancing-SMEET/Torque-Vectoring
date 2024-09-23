% Load data

Md = load("Md.mat");
time = Md.data.Time;
values = Md.data.Data;


% Plotting the data

plot(time, values);


% Plot the maximum difference

xlabel('t[s]', 'FontSize', 20);
ylabel('Md[Nm]', 'FontSize', 20);

% range
xlim([0,30])

set(gca, 'FontSize', 15); % Set font size for axis tick labels
set(gcf, 'Position', [100, 100, 1000, 500]);
% Add legend
legend('Md','FontSize', 15);
saveas(figure(1), 'C:\Users\jm538\Desktop\Code\TV_IITP\Paper\Figures\Md.tiff');
%% Sinus Test
% Data import
sinus_TV = load('C:\Users\jm538\Desktop\Code\TV_IITP\Simulink\results\TV\sinus\Sinus_test_TV.mat');
sinus_No = load('C:\Users\jm538\Desktop\Code\TV_IITP\Simulink\results\TV\sinus\Sinus_test_No.mat');

time = sinus_TV.data{1}.Values.Time;
time_No = sinus_No.data.Time;
r_des = sinus_TV.data{1}.Values.Data;
r_TV = sinus_TV.data{2}.Values.Data;
r_No = sinus_No.data.Data;

% plot
figure
plot(time,r_des,"LineWidth", 1);

hold on

plot(time,r_TV,"LineWidth", 1)

hold on

plot(time_No,r_No,"LineWidth", 1)
xlim([0,30])
fontSize = 18;
legendFontSize = 12;
set(gca,'FontSize',15)
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('yaw rate (rad/s)', 'FontSize', fontSize);
legend('{\gamma}_{desired}', '{\gamma}_{with TV}','{\gamma}_{without TV}', 'FontSize', legendFontSize);
set(gcf, 'Position', [100, 100, 1000, 500]);





%% Steady Test
% Data import
steady_TV = load('C:\Users\jm538\Desktop\Code\TV_IITP\Simulink\results\TV\steady\Steady_test_TV.mat');
steady_No = load('C:\Users\jm538\Desktop\Code\TV_IITP\Simulink\results\TV\steady\Steady_test.mat');

time = steady_TV.data{1}.Values.Time;
time_No = steady_No.data{2}.Values.Time;
r_des_TV = steady_TV.data{1}.Values.Data;
r_des_No = steady_No.data{2}.Values.Data;
r_TV = steady_TV.data{2}.Values.Data;
r_No = steady_No.data{1}.Values.Data;

% plot
figure
plot(time_No,r_des_No,"LineWidth", 1);

hold on
plot(time,r_TV,"LineWidth", 1)
hold on
plot(time_No,r_No,"LineWidth", 1)

fontSize = 18;
legendFontSize = 12;

set(gca,'FontSize',15)
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('yaw rate (rad/s)', 'FontSize', fontSize);
legend('{\gamma}_{desired}', '{\gamma}_{with TV}', '{\gamma}_{without TV}', 'FontSize', legendFontSize);
set(gcf, 'Position', [100, 100, 1000, 500]);

xlim([0,30])

figure
plot(time_No,r_des_No,"LineWidth", 1);

hold on
plot(time,r_TV,"LineWidth", 1)
hold on
plot(time_No,r_No,"LineWidth", 1)

fontSize = 18;
legendFontSize = 12;

set(gca,'FontSize',15)
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('yaw rate (rad/s)', 'FontSize', fontSize);
legend('{\gamma}_{desired}', '{\gamma}_{with TV}', '{\gamma}_{without TV}', 'FontSize', legendFontSize);
set(gcf, 'Position', [100, 100, 1000, 500]);

xlim([10,30])



%% race track

% Data import
race_TV = load('C:\Users\jm538\Desktop\Code\TV_IITP\Simulink\results\TV\race\race_tv.mat');
race_No = load('C:\Users\jm538\Desktop\Code\TV_IITP\Simulink\results\TV\race\race_No.mat');

time = race_TV.data{1}.Values.Time;
time_No = race_No.data.Time;
r_des_TV = race_TV.data{1}.Values.Data;
r_TV = race_TV.data{2}.Values.Data;
r_No = race_No.data.Data;

% plot

plot(time,r_des_TV);
hold on
plot(time,r_TV)
hold on

plot(time_No,r_No)


fontSize = 18;
legendFontSize = 12;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('raw rate (rad/s)', 'FontSize', fontSize);
legend('yaw rate_{desired}', 'yaw rate_{torque vectoring}','yaw rate_{base}', 'FontSize', legendFontSize);


