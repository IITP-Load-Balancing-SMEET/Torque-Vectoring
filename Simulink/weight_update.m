function [sys, x0, str, ts] = weight_update(t, x, u, flag)
    switch flag
        case 0
            [sys, x0, str, ts] = mdlInitializeSizes;
        case 2
            sys = mdlUpdate(t, x, u);
        case 3
            sys = mdlOutputs(t, x, u);
        case {1, 4, 9}
            sys = [];
        otherwise
            error(['Unhandled flag = ', num2str(flag)]);
    end
end

function [sys, x0, str, ts] = mdlInitializeSizes()
    % Initialization
    sizes = simsizes;
    sizes.NumContStates  = 0;
    sizes.NumDiscStates  = 3; % State 1: weight, State 2: previous value, State 3: switch
    sizes.NumOutputs     = 1; % Output the weight
    sizes.NumInputs      = 2; % Input: current value only
    sizes.DirFeedthrough = 0; % No direct feedthrough
    sizes.NumSampleTimes = 1; % One sample time

    sys = simsizes(sizes);
    
    % Initial conditions
    initialWeight = 1;        % Initial weight value
    initialPrevValue = 0;
    switch_value = 0; % Initial previous value
    
    x0 = [initialWeight; initialPrevValue; switch_value]; % Initial states
    str = []; % No state ordering
    ts = [0 0]; % Sample time: [period, offset]
end

function sys = mdlUpdate(~, x, u)
    alpha = 0.01;
    threshold = 0.05;


    % Get the current state values
    weight = x(1);
    prev_value = x(2);
    prev_switch = x(3);

    
    % Get the current input (current value)
    current_value = u(1);

    % Update the weight based on the comparison
    if prev_switch ~= u(2)
        weight = 1;
    end
    
    if current_value * prev_value < 0
        weight = weight * (1-alpha);
    elseif abs(current_value) < threshold
         weight = weight * (1-alpha); % Decrement weight
    else
        if abs(current_value) > abs(prev_value)
            weight = weight * (1+alpha); % Increment weight
            if weight > 2
                weight = 2;
            end
        end
    end
    
    % Update the states
    sys = [weight; current_value; u(2)]; % Updated weight and store current value as the previous value for next step
end

function sys = mdlOutputs(~, x, ~)
    % Output the current weight
    sys = x(1); % Output the weight (first state)
end
