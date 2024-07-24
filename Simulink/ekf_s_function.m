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

function [Fy_b, a, Fz] = dugoff(Vx, Vy, delta1, delta2, delta3, delta4, m, gamma, lf, lr, tr, h_cg, ax, ay)
    global delta1_prev delta2_prev delta3_prev delta4_prev
    global ax_prev ay_prev
    
    mu = 1.0; % road friction coefficient   
    load = m * 9.806; % Total vertical load
    L = lf + lr;
    ht =  0.5 * tr;

    if (isnan(delta1) || isnan(delta2) || isnan(delta3) || isnan(delta4) || isnan(ax_prev) || isnan(ay_prev))
        delta1 = delta1_prev;
        delta2 = delta2_prev;
        delta3 = delta3_prev;
        delta4 = delta4_prev;
        ax = ax_prev;
        ay = ay_prev;
    end

    if abs(Vx) <= 0.25 % prevent divided by zero (standstill velocity)
        Fy_b = [0.0; 
                0.0; 
                0.0; 
                0.0];

        a = [0.0; 
             0.0; 
             0.0; 
             0.0];

        Fz = [load*(lr / (2*L)) - (m*ax*h_cg / (2*L)) + (m*ay*h_cg*lr / (ht*L));
              load*(lr / (2*L)) - (m*ax*h_cg / (2*L)) - (m*ay*h_cg*lr / (ht*L));
              load*(lf / (2*L)) + (m*ax*h_cg / (2*L)) + (m*ay*h_cg*lf / (ht*L));
              load*(lf / (2*L)) + (m*ax*h_cg / (2*L)) - (m*ay*h_cg*lf / (ht*L))];

        delta1_prev = delta1;
        delta2_prev = delta2;
        delta3_prev = delta3;
        delta4_prev = delta4;
        ax_prev = ax;
        ay_prev = ay;
        
        return;
    
    else
        Fz = [load*(lr / (2*L)) - (m*ax*h_cg / (2*L)) + (m*ay*h_cg*lr / (ht*L));
              load*(lr / (2*L)) - (m*ax*h_cg / (2*L)) - (m*ay*h_cg*lr / (ht*L));
              load*(lf / (2*L)) + (m*ax*h_cg / (2*L)) + (m*ay*h_cg*lf / (ht*L));
              load*(lf / (2*L)) + (m*ax*h_cg / (2*L)) - (m*ay*h_cg*lf / (ht*L))];   

        a = [-delta1 + atan((Vy + lf*gamma) / (Vx)); 
             -delta2 + atan((Vy + lf*gamma) / (Vx)); 
             -delta3 + atan((Vy - lr*gamma) / (Vx));
             -delta4 + atan((Vy - lr*gamma) / (Vx))];
        
        Cy = (-0.005789 .* Fz.^2)  + (3.500753 .* Fz) + 30000;
        lambda = (mu * Fz) ./ (2.0 * Cy .* abs(tan(a)));
        flambda = ones(size(lambda));
        mask = lambda < 1.0;
        flambda(mask) = (2 - lambda(mask)) .* lambda(mask); 
    
        Fy_b = [-Cy(1) * tan(a(1)) * flambda(1);
                -Cy(2) * tan(a(2)) * flambda(2);
                -Cy(3) * tan(a(3)) * flambda(3);
                -Cy(4) * tan(a(4)) * flambda(4)];

        delta1_prev = delta1;
        delta2_prev = delta2;
        delta3_prev = delta3;
        delta4_prev = delta4;
        ax_prev = ax;
        ay_prev = ay;
    end
end

function [x_pred, Fy_b] = stateTransitionFunction(x_est, u_c, params, ax, ay, dt)
    Vx = x_est(1);
    Vy = x_est(2);
    gamma = x_est(3);

    Fy_fl = x_est(4);
    Fy_fr = x_est(5);
    Fy_rl = x_est(6);
    Fy_rr = x_est(7);
    
    delta1 = u_c(1);
    delta2 = u_c(2);
    delta3 = u_c(3);
    delta4 = u_c(4);
    delta = (delta1 + delta2) / 2;
    Fx_fl = u_c(5);
    Fx_fr = u_c(6);
    Fx_rl = u_c(7); 
    Fx_rr = u_c(8);

    m = params(1);
    Iz = params(2);
    radius = params(3);
    Cav = params(4);
    lf = params(5);
    lr = params(6);
    tr = params(7);
    sigma = params(8);
    h_cg = params(9);

    [Fy_b, ~] = dugoff(Vx, Vy, delta1, delta2, delta3, delta4, m, gamma, lf, lr, tr, h_cg, ax, ay);

    x1_dot = 1/m * (((Fx_fl + Fx_fr)*cos(delta)) - ((Fy_fl + Fy_fr)*sin(delta)) + Fx_rl + Fx_rr - (Cav*Vx^2)) + (Vy*gamma); 
    x2_dot = 1/m * (((Fx_fl + Fx_fr)*sin(delta)) + ((Fy_fl + Fy_fr)*cos(delta)) + Fy_rl + Fy_rr) - (Vx*gamma);
    x3_dot = 1/Iz * (lf * (((Fx_fl + Fx_fr)*sin(delta)) + ((Fy_fl + Fy_fr)*cos(delta))) + ...
        (tr * (((Fx_fl - Fx_fr)*cos(delta)) + ((-Fy_fl + Fy_fr)*sin(delta)))) - ...
        (lr * (Fy_rl + Fy_rr)));

    x4_dot = (Vx/sigma) * (-Fy_fl + Fy_b(1)); 
    x5_dot = (Vx/sigma) * (-Fy_fr + Fy_b(2));  
    x6_dot = (Vx/sigma) * (-Fy_rl + Fy_b(3));  
    x7_dot = (Vx/sigma) * (-Fy_rr + Fy_b(4)); %  have to check what kind of values to insert  
    x_pred = x_est + [x1_dot; x2_dot; x3_dot; x4_dot; x5_dot; x6_dot; x7_dot] * dt;  

end

function [sys, x0, str, ts] = mdlInitializeSizes
    global Q R
    global delta1_prev delta2_prev delta3_prev delta4_prev
    global ax_prev ay_prev
    global t_prev

    sizes = simsizes;
    sizes.NumContStates  = 0; % 연속 상태의 수
    sizes.NumDiscStates  = 7 * 7 + 7;  % 개별상태의 수, [x; P(:)] where P is 7x7, thus 49 elements
    sizes.NumOutputs     = 15;  % 출력 수, State estimate [V_x; V_y; gamma; F_yfl; F_yfr; F_yrl; F_yrr] and vector 'a', "Fz"
    sizes.NumInputs      = 13;  % 입력 수, [delta1; delta2; T_d1 - T_b1; T_d2 - T_b2; T_d3 - T_b3; T_d4 - T_b4; Z (Vx, Vy, yaw_rate, ax, ay)]
    sizes.DirFeedthrough = 0; % reserved
    sizes.NumSampleTimes = 1; % 샘플 횟수
    
    sys = simsizes(sizes);

    P0 = diag([1, 1, 0.001, 10, 10, 10, 10]);
    x0  = [0.0; 0.0; 0.0; 0.0; 0.0; 0.0; 0.0; reshape(P0, [], 1)];  % Initial state and covariance matrix
    str = [];
    ts  = [0.001 0];  % Sample time
    
    Q = diag([1, 1, 0.001, 1, 1, 1, 1]); % Process Noise Covariance Matrix
    R = diag([0.0001, 0.0001, 0.000001, 0.0001, 0.0001]); % Measurement Noise Covariance Matrix 

    delta1_prev = 0.0;
    delta2_prev = 0.0;
    delta3_prev = 0.0;
    delta4_prev = 0.0;
    ax_prev = 0.0;
    ay_prev = 0.0;
    t_prev = 0.0;
end

function sys = mdlDerivatives(t, x, u)
    sys = [];
end

function sys = mdlUpdate(t, x, u)
    global Q R
    global t_prev

    % Define system parameters
    m = (226.26 + 70);  % Mass of the vehicle + 70 kg Load
    Iz = 146.827; % moment of inertia about yaw axis
    radius = 0.262;  % Wheel radius
    Cav = (0.5 * 1.205 * 3.237 * 1.1);  % Aerodynamic coefficient (0.5 * air_density(=1.205 kg m−3) * air_drag * cross_sectional_area)
    lf = 0.813;  % Distance from CG to front axle
    lr = 0.787;  % Distance from CG to rear axle
    tr = 1.242;  % Track width
    sigma = 0.1;
    h_cg = 0.3; % height from ground to C.G

    params = [m, Iz, radius, Cav, lf, lr, tr, sigma, h_cg];
    
    dt = (t - t_prev);

    % Unpack state vector
    x_est = x(1:7);
    P_est = reshape(x(8:end), 7, 7);

    u_c = u(1:8); % control input
    delta = (u_c(1) + u_c(2)) / 2; % average steer angle
    z = u(9:end); % measurement

    ax = z(4);
    ay = z(5);
    
    [x_pred, Fy_b] = stateTransitionFunction(x_est, u_c, params, ax, ay, dt);

    % Jacobian
    F = [-(2*Cav*x_est(1))/m, x_est(3), x_est(2), -sin(delta)/m, -sin(delta)/m, 0, 0;
         -x_est(3), 0, -x_est(1), cos(delta)/m, cos(delta)/m, 1/m, 1/m;
         0, 0, 0, (lf*cos(delta) - tr*sin(delta))/Iz, (lf*cos(delta) + tr*sin(delta))/Iz, -lr/Iz, -lr/Iz;
         (Fy_b(1) - x_est(4))/sigma, 0, 0, -x_est(1)/sigma, 0, 0, 0;
         (Fy_b(2) - x_est(5))/sigma, 0, 0, 0, -x_est(1)/sigma, 0, 0;
         (Fy_b(3) - x_est(6))/sigma, 0, 0, 0, 0, -x_est(1)/sigma, 0;
         (Fy_b(4) - x_est(7))/sigma, 0, 0, 0, 0, 0, -x_est(1)/sigma] * dt;

    A = eye(7) + F;
    
    P_pred = (A*P_est*A') + Q;

    % Kalman Gain
    h_prior = [x_pred(1);
               x_pred(2);
               x_pred(3);
               (1/m) * (((u(5) + u(6))*cos(delta)) - ((x_pred(4) + x_pred(5))*sin(delta)) + (u(7) + u(8)) - (Cav*x_pred(1)^2));
               (1/m) * (((u(5) + u(6))*sin(delta)) + ((x_pred(4) + x_pred(5))*cos(delta)) + x_pred(6) + x_pred(7))];
    
    H = [1, 0, 0, 0, 0, 0, 0;
         0, 1, 0, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0, 0;
         -2*Cav*x_est(1)/m, 0, 0, -sin(delta)/m, -sin(delta)/m, 0, 0;
         0, 0, 0, cos(delta)/m, cos(delta)/m, 1/m, 1/m];

    K = P_pred * H' / (H*P_pred*H' + R);

    % Update step
    alpha = 0.75;
    innov = z - h_prior;
    Q = (alpha * Q) + ((1-alpha) * (K*innov*innov'*K'));

    x_est = x_pred + (K * innov);
    P_est = (eye(7) - K * H) * P_pred;
    
    sys = [x_est; P_est(:)];
    t_prev = t;
end

function sys = mdlOutputs(t, x, u)
    x_est = x(1:7);
    
    % Calculate vector 'a' based on the current state
    Vx = x_est(1);
    Vy = x_est(2);
    gamma = x_est(3);

    delta1 = u(1);
    delta2 = u(2);
    delta3 = u(3);
    delta4 = u(4);

    m = 296.26;  % Update as per your specific parameters
    lf = 0.813;
    lr = 0.787;
    tr = 1.242;
    h_cg = 0.3;

    ax = u(10);
    ay = u(11);

    [~, a, Fz] = dugoff(Vx, Vy, delta1, delta2, delta3, delta4, m, gamma, lf, lr, tr, h_cg, ax, ay);
    % Output the state estimate and the vector 'a'
    sys = [x_est; a; Fz];
end
