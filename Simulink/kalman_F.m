function [sys, x0, str, ts] = ekf_sfunc(t, x, u, flag, f_func, F_jacob, h_func, H_jacob, Q, R, dt)
    switch flag
        case 0
            [sys, x0, str, ts] = mdlInitializeSizes;
        case 1
            sys = mdlDerivatives(t, x, u);
        case 3
            sys = mdlOutputs(t, x, u, f_func, F_jacob, h_func, H_jacob, Q, R, dt);
        case {2, 4, 9}
            sys = [];
        otherwise
            error(['Unhandled flag = ', num2str(flag)]);
    end
end

function [sys, x0, str, ts] = mdlInitializeSizes
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 0;
    sizes.NumOutputs     = 2; % [x_upd; P_upd] concatenated
    sizes.NumInputs      = 3; % [x; P; z; u]
    sizes.DirFeedthrough = 1;
    sizes.NumSampleTimes = 1;
    sys = simsizes(sizes);
    x0  = [];
    str = [];
    ts  = [dt 0];
end

function sys = mdlDerivatives(~, ~, ~)
    sys = [];
end

function sys = mdlOutputs(~, ~, u, f_func, F_jacob, h_func, H_jacob, Q, R, dt)
    % Extract inputs
    x = u(1:end/3);
    P = reshape(u(end/3+1:2*end/3), sqrt(numel(u(end/3+1:2*end/3))), []);
    z = u(2*end/3+1:end);
    control_input = u(end);

    % EKF Prediction
    [x_pred, P_pred] = ekf_predict(x, P, control_input, Q, f_func, F_jacob, dt);

    % EKF Update
    [x_upd, P_upd] = ekf_update(x_pred, P_pred, z, R, h_func, H_jacob);

    % Output updated state and covariance
    sys = [x_upd; P_upd(:)];
end
