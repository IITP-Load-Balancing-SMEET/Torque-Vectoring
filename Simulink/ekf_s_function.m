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
    sizes.NumInputs      = 11;  % 입력 수, [delta1; delta2; T_d1 - T_b1; T_d2 - T_b2; T_d3 - T_b3; T_d4 - T_b4; Z (Vx, Vy, yaw_rate, ax, ay)]
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

function [a_fl, a_fr, a_rl, a_rr] = tire_slip(Vx, Vy, gamma, lf, lr, delta1, delta2)
    if Vx < 1.0 % prevent divided by zero
        a_fl = 0.0;
        a_fr = 0.0;
        a_rl = 0.0;
        a_rr = 0.0;
    else
        a_fl = -delta1 + atan2((Vy + lf*gamma), Vx);
        a_fr = -delta2 + atan2((Vy + lf*gamma), Vx);
        a_rl = atan2((Vy - lr*gamma), Vx);
        a_rr = a_rl;
    end
end

function [Cy_fl, Cy_fr, Cy_rl, Cy_rr] = cornering_stiffness(Fy_fl, Fy_fr, Fy_rl, Fy_rr)
    if (a_fl == 0.0) || (a_fr == 0.0) || (a_rl == 0.0) || (a_rr == 0.0)
        Cy_fl = 50000; % need to be fixed
        Cy_fr = 50000; % need to be fixed
        Cy_rl = 50000; % need to be fixed
        Cy_rr = 50000; % need to be fixed
    else
        Cy_fl = Fy_fl / a_fl;
        Cy_fr = Fy_fr / a_fr;
        Cy_rl = Fy_rl / a_rl;
        Cy_rr = Fy_rr / a_rr;
    end
end

function [Fy_flb, Fy_frb, Fy_rlb, Fy_rrb] = dugoff(m, Vx, Vy, gamma, lf, lr, delta1, delta2, Fy_fl, Fy_fr, Fy_rl, Fy_rr)
    [a_fl, a_fr, a_rl, a_rr] = tire_slip(Vx, Vy, gamma, lf, lr, delta1, delta2);
    [Cy_fl, Cy_fr, Cy_rl, Cy_rr] = cornering_stiffness(Fy_fl, Fy_fr, Fy_rl, Fy_rr);
    
    Fz = m * -9.81; % Total vertical load
    mu = 0.9; % road firction coefficient

    lamda_fl = (mu * Fz) / (2 * Cy_fl * abs(tan(a_fl)));
    lamda_fr = (mu * Fz) / (2 * Cy_fr * abs(tan(a_fr)));
    lamda_rl = (mu * Fz) / (2 * Cy_rl * abs(tan(a_rl)));
    lamda_rr = (mu * Fz) / (2 * Cy_rr * abs(tan(a_rr)));
    
    lambda = [lambda_fl, lambda_fr, lambda_rl, lambda_rr];

    for i = 1:length(lambda)
        if lambda(i) >= 1
            lambda(i) = 1.0;
        else
            lambda(i) = (2 - lambda(i)) * lambda(i);
        end    
    end

    Fy_flb = -Cy_fl * tan(a_fl) * lambda(1);
    Fy_frb = -Cy_fr * tan(a_fr) * lambda(2);
    Fy_rlb = -Cy_rl * tan(a_rl) * lambda(3);
    Fy_rrb = -Cy_rr * tan(a_rr) * lambda(4);
end



function sys = mdlUpdate(t, x, u)
    % Unpack state vector
    x_est = x(1:7);
    P_est = reshape(x(8:end), 7, 7);

    % Control input and measurement
    delta1 = u(1); % left steer
    delta2 = u(2); % right steer
    delta = mean(u(1:2)); % average steer
    torque = u(3:6);
    z = u(7:end);

    % Define system parameters
    dt = 0.001;
    m = (226.26 + 70);  % Mass of the vehicle + 70 kg Load
    radius = 0.262;  % Wheel radius
    Cav = (0.5 * 1.293 * 1.4285 * 1.1);  % Aerodynamic coefficient (0.5 * air_density(=1.293 kg m−3) * air_drag * cross_sectional_area)
    lf = 0.813;  % Distance from CG to front axle
    lr = 0.787;  % Distance from CG to rear axle
    track = 1.242;  % Track width
    sigma = 0.1;  % relaxation Length to tire force (need to be fixed)

    Q = diag([1, 1, 0.0001, 1000, 1000, 1000, 1000]); % Process Noise Covariance Matrix
    R = diag([0.0001, 0.001, 0.000001, 0.0001, 0.0001]); % Measurement Noise Covariance Matrix 
    
    [x_pred, Fy_b] = stateTransitionFunction(x_est, delta, delta1, delta2, torque, m, radius, Cav, lf, lr, track, dt, sigma, F_yfl_max);
    
    % Jacobian
    F = [-(2*Cav*x_est(1))/m, x_est(3), x_est(2), -sin(delta)/m, -sin(delta)/m, 0, 0;
         x_est(3), 0, x_est(1), cos(delta)/m, cos(delta)/m, 1/m, 1/m;
         0, 0, 0, (lf*cos(delta) - track*cos(delta))/m, (lf*cos(delta) + track*cos(delta))/m, -lr/m, -lr/m;
         (Fy_b(1) - x_est(4))/sigma, 0, 0, -x_est(1)/sigma, 0, 0, 0;
         (Fy_b(2) - x_est(5))/sigma, 0, 0, 0, -x_est(1)/sigma, 0, 0;
         (Fy_b(3) - x_est(6))/sigma, 0, 0, 0, 0, -x_est(1)/sigma, 0;
         (Fy_b(4) - x_est(7))/sigma, 0, 0, 0, 0, 0, -x_est(1)/sigma];
    
    P_pred = F*P_est*F' + Q;

    % Kalman Gain
    H = [1, 0, 0, 0, 0, 0, 0;
         0, 1, 0, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0, 0;
         -2*Cav*x_est(1)/m, 0, 0, -sin(delta)/m, -sin(delta)/m, 0, 0;
         0, 0, 0, cos(delta)/m, cos(delta)/m, 1/m, 1/m];

    K = P_pred * H' / (H*P_pred*H' + R);
    h = [x_pred(1);
         x_pred(2);
         x_pred(3);
         (1/m) * ( (u(2) + u(3))/radius * cos(delta) - (x_pred(4) + x_pred(5)) * sin(delta) + (u(4) + u(5))/radius - Cav * x_pred(1)^2 );
         (1/m) * ( (u(2) + u(3))/radius * sin(delta) + (x_pred(4) + x_pred(5)) * cos(delta) + x_pred(6) + x_pred(7))];

    % Update step
    x_est = x_pred + (K * (z - h));
    P_est = (eye(7) - K * H) * P_pred;

    % Pack state vector
    sys = [x_est; P_est(:)];
end

function sys = mdlOutputs(t, x, u)
    sys = x(1:7);  % Output the state estimate
end

function [x_pred, Fy_b] = stateTransitionFunction(x, delta, delta1, delta2, torque, m, radius, C_av, l_f, l_r, track, dt, sigma)
    % Define the state transition function f(x)
    Vx = x(1);
    Vy = x(2);
    gamma = x(3);
    Fy_fl = x(4);
    Fy_fr = x(5);
    Fy_rl = x(6);
    Fy_rr = x(7);

    T_fl = torque(1);
    T_fr = torque(2);
    T_rl = torque(3); 
    T_rr = torque(4);
    
    [Fy_flb, Fy_frb, Fy_rlb, Fy_rrb] = dugoff(m, Vx, Vy, gamma, lf, lr, delta1, delta2, Fy_fl, Fy_fr, Fy_rl, Fy_rr);

    x1_dot = 1/m * (1/radius * (T_fl + T_fr) * cos(delta) - (Fy_fl + Fy_fr) * sin(delta) + 1/radius * (T_rl + T_rr) - C_av * Vx^2) + Vy * gamma;
    x2_dot = 1/m * (1/radius * (T_fl + T_fr) * sin(delta) + (Fy_fl + Fy_fr) * cos(delta) + (Fy_rl + Fy_rr)) + Vx * gamma;
    x3_dot = 1/m * (l_f * (1/radius * (T_fl + T_fr) * sin(delta) + (Fy_fl + Fy_fr) * cos(delta)) +  track * (1/radius * (T_fl - T_fr) * cos(delta) + (-Fy_fl + Fy_fr) * cos(delta) + 1/radius * (T_rl - T_rr)) - l_r * (Fy_rl + Fy_rr));
    x4_dot = Vx / sigma * (-Fy_fl + Fy_flb);
    x5_dot = Vx / sigma * (-Fy_fr + Fy_frb);  
    x6_dot = Vx / sigma * (-Fy_rl + Fy_rlb);  
    x7_dot = Vx / sigma * (-Fy_rr + Fy_rrb); %  have to check what kind of values to insert  
 
    x_pred = x + [x1_dot; x2_dot; x3_dot; x4_dot; x5_dot; x6_dot; x7_dot] * dt;  
    Fy_b = [Fy_flb, Fy_frb, Fy_rlb, Fy_rrb];
end
