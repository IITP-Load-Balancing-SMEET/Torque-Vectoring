function [sys, x0, str, ts] = ekf_s_function(t, x, u, flag)
    switch flag
        case 0
            [sys, x0, str, ts] = mdlInitializeSizes;
        case 1
            sys = mdlDerivatives(t, x, u);
        case 2
            sys = mdlUpdate(t, x, u);
        case 3
            sys = mdlOutputs(t, x, u);
        case {4, 9}
            sys = [];
        otherwise
            error(['Unhandled flag = ', num2str(flag)]);
    end
end

function [sys, x0, str, ts] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0; % 연속 상태의 수
    sizes.NumDiscStates  = 7 * 7 + 7;  % 개별상태의 수, [x; P(:)] where P is 7x7, thus 49 elements
    sizes.NumOutputs     = 7;  % 출력 수, State estimate [V_x; V_y; gamma; F_yfl; F_yfr; F_yrl; F_yrr]
    sizes.NumInputs      = 10;  % 입력 수, [delta; T_d1 - T_b1; T_d2 - T_b2; T_d3 - T_b3; T_d4 - T_b4; Z (Vx, Vy, yaw_rate, ax, ay)]
    sizes.DirFeedthrough = 0; % reserved
    sizes.NumSampleTimes = 1; % 샘플 횟수
    
    sys = simsizes(sizes);

    P0 = diag([1, 1, 0.0001, 1000, 1000, 1000, 1000]);
    x0  = [0; 0; 0; 0; 0; 0; 0; reshape(P0, [], 1)];  % Initial state and covariance matrix
    str = [];
    ts  = [-1 0];  % Sample time
end

function sys = mdlDerivatives(t, x, u)
    sys = [];
end

function sys = mdlUpdate(t, x, u)
    % Unpack state vector
    x_est = x(1:7);
    P_est = reshape(x(8:end), 7, 7);

    % Control input and measurement
    delta = u(1);
    Torque = u(2:5);
    z = u(6:end);

    % Define system parameters
    m = (226.26 + 70);  % Mass of the vehicle + 70 kg Load
    Radius = 0.262;  % Wheel radius
    C_av = (0.5 * 1.293 * 1.4285 * 1.1);  % Aerodynamic coefficient (0.5 * air_density(=1.293 kg m−3) * air_drag * cross_sectional_area)
    l_f = 0.813;  % Distance from CG to front axle
    l_r = 0.787;  % Distance from CG to rear axle
    track = 1.242;  % Track width
    sigma = 0.1;  % relaxation Length to tire force
    F_yfl_max = 1000;  % Maximum lateral tire force

    Q = diag([1, 1, 0.0001, 1000, 1000, 1000, 1000]); % Process Noise Covariance Matrix
    R = diag([0.0001, 0.001, 0.000001, 0.0001, 0.0001]); % Measurement Noise Covariance Matrix 

    % Jacobian
    F = [-(2*C_av*x_est(1))/m, x_est(3), x_est(2), -sin(delta)/m, -sin(delta)/m, 0, 0;
         x_est(3), 0, x_est(1), cos(delta)/m, cos(delta)/m, 1/m, 1/m;
         0, 0, 0, (l_f*cos(delta) - track*cos(delta))/m, (l_f*cos(delta) + track*cos(delta))/m, -l_r/m, -l_r/m;
         (F_yfl_max - x_est(4))/sigma, 0, 0, -x_est(1)/sigma, 0, 0, 0;
         (F_yfl_max - x_est(5))/sigma, 0, 0, 0, -x_est(1)/sigma, 0, 0;
         (F_yfl_max - x_est(6))/sigma, 0, 0, 0, 0, -x_est(1)/sigma, 0;
         (F_yfl_max - x_est(7))/sigma, 0, 0, 0, 0, 0, -x_est(1)/sigma];

    H = [1, 0, 0, 0, 0, 0, 0;
         0, 1, 0, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0, 0;
         -2*C_av*x_est(1)/m, 0, 0, -sin(delta)/m, -sin(delta)/m, 0, 0;
         0, 0, 0, cos(delta)/m, cos(delta)/m, 1/m, 1/m];
    
    dt = 0.001;
    
    x_pred = stateTransitionFunction(x_est, delta, Torque, m, Radius, C_av, l_f, l_r, track, dt, sigma, F_yfl_max);
    P_pred = F*P_est*F' + Q;

    % Kalman Gain
    K = P_pred * H' / (H*P_pred*H' + R);

    h = [x_pred(1);
         x_pred(2);
         x_pred(3);
         (1/m) * ( (u(2) + u(3))/Radius * cos(delta) - (x_pred(4) + x_pred(5)) * sin(delta) + (u(4) + u(5))/Radius - C_av * x_pred(1)^2 );
         (1/m) * ( (u(2) + u(3))/Radius * sin(delta) + (x_pred(4) + x_pred(5)) * cos(delta) + x_pred(6) + x_pred(7))];

    % Update step
    %x_est = x_pred + (K * (z - h));
    %P_est = (eye(7) - K * H) * P_pred;
    x_est = x_pred;
    P_est = P_pred;

    disp(x_est);

    % Pack state vector
    sys = [x_est; P_est(:)];
end

function sys = mdlOutputs(t, x, u)
    sys = x(1:7);  % Output the state estimate
end

function x_pred = stateTransitionFunction(x, delta, Torque, m, R, C_av, l_f, l_r, track, dt, sigma, F_yfl_max)
    % Define the state transition function f(x)
    V_x = x(1);
    V_y = x(2);
    gamma = x(3);
    F_yfl = x(4);
    F_yfr = x(5);
    F_yrl = x(6);
    F_yrr = x(7);

    T_fl = Torque(1);
    T_fr = Torque(2);
    T_rl = Torque(3); 
    T_rr = Torque(4);

    x1_dot = 1/m * (1/R * (T_fl + T_fr) * cos(delta) - (F_yfl + F_yfr) * sin(delta) + 1/R * (T_rl + T_rr) - C_av * V_x^2) + V_y * gamma;
    x2_dot = 1/m * (1/R * (T_fl + T_fr) * sin(delta) + (F_yfl + F_yfr) * cos(delta) + (F_yrl + F_yrr)) + V_x * gamma;
    x3_dot = 1/m * (l_f * (1/R * (T_fl + T_fr) * sin(delta) + (F_yfl + F_yfr) * cos(delta)) +  track * (1/R * (T_fl - T_fr) * cos(delta) + (-F_yfl + F_yfr) * cos(delta) + 1/R * (T_rl - T_rr)) - l_r * (F_yrl + F_yrr));
    x4_dot = V_x / sigma * (-F_yfl + F_yfl_max);
    x5_dot = V_x / sigma * (-F_yfr + F_yfl_max);  
    x6_dot = V_x / sigma * (-F_yrl + F_yfl_max);  
    x7_dot = V_x / sigma * (-F_yrr + F_yfl_max); %  have to check what kind of values to insert  
 
    x_pred = x + [x1_dot; x2_dot; x3_dot; x4_dot; x5_dot; x6_dot; x7_dot] * dt;  
end
