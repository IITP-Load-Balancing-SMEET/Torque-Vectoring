% Load data

Md = load("Md.mat");
time = Md.data.Time;
values = Md.data.Data;

time = time(50003:150000);
values = values(50003:150000);

time = time - time(1);
% Plotting the data

plot(time, values);


% Plot the maximum difference

xlabel('t[s]', 'FontSize', 20);
ylabel('Md[Nm]', 'FontSize', 20);

% range
xlim([0,100])

set(gca, 'FontSize', 15); % Set font size for axis tick labels

% Add legend
legend('Md','FontSize', 15);

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

plot(time,r_des);

hold on

plot(time,r_TV)

fontSize = 18;
legendFontSize = 12;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('raw rate (rad/s)', 'FontSize', fontSize);
legend('yaw rate_{desired}', 'yaw rate_{torque vectoring}', 'FontSize', legendFontSize);
figure
plot(time,r_des);
hold on
plot(time_No,r_No)
fontSize = 18;
legendFontSize = 12;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('raw rate (rad/s)', 'FontSize', fontSize);
legend('yaw rate_{desired}','yaw rate_{base}', 'FontSize', legendFontSize);



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

plot(time,r_des_TV);

hold on
plot(time,r_TV)

fontSize = 18;
legendFontSize = 12;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('raw rate (rad/s)', 'FontSize', fontSize);
legend('yaw rate_{desired}', 'yaw rate_{torque vectoring}', 'FontSize', legendFontSize);


plot(time_No,r_des_No)
hold on
plot(time_No,r_No)

fontSize = 18;
legendFontSize = 12;
xlabel('Time (s)', 'FontSize', fontSize);
ylabel('raw rate (rad/s)', 'FontSize', fontSize);
legend('yaw rate_{desired}', 'yaw rate_{base}', 'FontSize', legendFontSize);



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


