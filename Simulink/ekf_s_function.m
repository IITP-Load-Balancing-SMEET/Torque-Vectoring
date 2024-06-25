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


function [Fy_flb, Fy_frb, Fy_rlb, Fy_rrb] = dugoff(Vx, Vy, delta1, delta2, m, gamma, lf, lr)
    % tire slip and cornering stiffness calculation

    if abs(Vx) <= 0.25 % prevent divided by zero
        Fy_flb = 0.01;
        Fy_frb = 0.01;
        Fy_rlb = 0.01;
        Fy_rrb = 0.01;

    else
        Fz = m * 9.81; % Total vertical load
        mu = 0.9; % road firction coefficient
        Cy = 3600;
        
        a_fl = -delta1 + atan((Vy + lf*gamma) / Vx);
        a_fr = -delta2 + atan((Vy + lf*gamma) / Vx);
        a_rl = atan((Vy - lr*gamma) / Vx);
        a_rr = atan((Vy - lr*gamma) / Vx);
            
        a = [a_fl, a_fr, a_rl, a_rr];

        lambda = (mu * Fz) ./ (2 * Cy .* abs(tan(a)));
        
        flambda = [0, 0, 0, 0];

        for i=1:length(lambda)
            if lambda(i) < 1.0
                flambda(i) = (2 - lambda(i)) * lambda(i);
            else
                flambda(i) = 1.0;
            end
        end
        
        Fy_flb = -Cy * tan(a_fl) * flambda(1);
        Fy_frb = -Cy * tan(a_fr) * flambda(2);
        Fy_rlb = -Cy * tan(a_rl) * flambda(3);
        Fy_rrb = -Cy * tan(a_rr) * flambda(4);
    end

end


function [x_pred, Fy_b] = stateTransitionFunction(x_est, u_c, params, dt)
    Vx = x_est(1);
    Vy = x_est(2);
    gamma = x_est(3);
    Fy_fl = x_est(4);
    Fy_fr = x_est(5);
    Fy_rl = x_est(6);
    Fy_rr = x_est(7);

    delta1 = u_c(1);
    delta2 = u_c(2);
    T_fl = u_c(3);
    T_fr = u_c(4);
    T_rl = u_c(5); 
    T_rr = u_c(6);

    m = params(1);
    Iz = params(2);
    radius = params(3);
    Cav = params(4);
    lf = params(5);
    lr = params(6);
    tr = params(7);
    sigma = params(8);
    
    Fx_fl = T_fl/radius; Fx_fr = T_fr/radius; Fx_rl = T_rl/radius; Fx_rr = T_rr/radius;

    [Fy_flb, Fy_frb, Fy_rlb, Fy_rrb] = dugoff(Vx, Vy, delta1, delta2, m, gamma, lf, lr);
    Fy_b = [Fy_flb, Fy_frb, Fy_rlb, Fy_rrb];

    x1_dot = 1/m * ((Fx_fl*cos(delta1)) + (Fx_fr*cos(delta2)) - (Fy_fl*sin(delta1)) - (Fy_fr*sin(delta2)) + Fx_rl + Fx_rr - (Cav*Vx^2)) + ...
        (Vy*gamma);
    x2_dot = 1/m * ((Fx_fl*sin(delta1)) + (Fx_fr*sin(delta2)) + (Fy_fl*cos(delta1)) + (Fy_fr*cos(delta2)) + Fy_rl + Fy_rr) - ...
        (Vx*gamma);
    x3_dot = lf/Iz * ((Fx_fl*sin(delta1) + Fx_fr*sin(delta2) + Fy_fl*cos(delta1) + Fy_fr*cos(delta2)) + ...
        (tr * (Fx_fl*cos(delta1) - Fx_fr*cos(delta2) - Fy_fl*sin(delta1) + Fy_fr*sin(delta2))) - ...
        (lr * (Fy_rl + Fy_rr)));
    x4_dot = (Vx/sigma) * (-Fy_fl + Fy_flb);
    x5_dot = (Vx/sigma) * (-Fy_fr + Fy_frb);  
    x6_dot = (Vx/sigma) * (-Fy_rl + Fy_rlb);  
    x7_dot = (Vx/sigma) * (-Fy_rr + Fy_rrb); %  have to check what kind of values to insert  

    x_pred = x_est + [x1_dot; x2_dot; x3_dot; x4_dot; x5_dot; x6_dot; x7_dot] * dt;  
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
    x0  = [0; 0; 0; 0.01; 0.01; 0.01; 0.01; reshape(P0, [], 1)];  % Initial state and covariance matrix
    str = [];
    ts  = [-1 0];  % Sample time
end


function sys = mdlDerivatives(t, x, u)
    sys = [];
end


function sys = mdlUpdate(t, x, u)
    % Define system parameters
    m = (226.26 + 70);  % Mass of the vehicle + 70 kg Load
    Iz = 146.827; % moment of inertia about yaw axis
    radius = 0.262;  % Wheel radius
    Cav = (0.5 * 1.293 * 3.237 * 1.1);  % Aerodynamic coefficient (0.5 * air_density(=1.293 kg m−3) * air_drag * cross_sectional_area)
    lf = 0.813;  % Distance from CG to front axle
    lr = 0.787;  % Distance from CG to rear axle
    tr = 1.242;  % Track width
    sigma = 0.1;  % relaxation Length to tire force (need to be fixed)
    params = [m, Iz, radius, Cav, lf, lr, tr, sigma];
    
    dt = 0.001; 

    % Unpack state vector
    x_est = x(1:7);
    P_est = reshape(x(8:end), 7, 7);

    u_c = u(1:6); % control input
    z = u(7:end); % measurement

    Q = diag([1, 1, 0.0001, 1000, 1000, 1000, 1000]); % Process Noise Covariance Matrix
    R = diag([0.0001, 0.001, 0.000001, 0.0001, 0.0001]); % Measurement Noise Covariance Matrix 
    
    [x_pred, Fy_b] = stateTransitionFunction(x_est, u_c, params, dt);
    

    delta = (u(1) + u(2)) / 2;

    % Jacobian
    F = [-(2*Cav*x_est(1))/m, x_est(3), x_est(2), -sin(delta)/m, -sin(delta)/m, 0, 0;
         -x_est(3), 0, -x_est(1), cos(delta)/m, cos(delta)/m, 1/m, 1/m;
         0, 0, 0, (lf*cos(delta) - tr*sin(delta))/Iz, (lf*cos(delta) + tr*sin(delta))/Iz, -lr/Iz, -lr/Iz;
         (Fy_b(1) - x_est(4))/sigma, 0, 0, -x_est(1)/sigma, 0, 0, 0;
         (Fy_b(2) - x_est(5))/sigma, 0, 0, 0, -x_est(1)/sigma, 0, 0;
         (Fy_b(3) - x_est(6))/sigma, 0, 0, 0, 0, -x_est(1)/sigma, 0;
         (Fy_b(4) - x_est(7))/sigma, 0, 0, 0, 0, 0, -x_est(1)/sigma] * dt;
    
    P_pred = (F*P_est*F') + Q;

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

    % 
    % ack state vector
    %x_est = x_pred;
    %P_est = P_pred;
    sys = [x_est; P_est(:)];
end


function sys = mdlOutputs(t, x, u)
    sys = x(1:7);  % Output the state estimate
end
