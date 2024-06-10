% Define the symbolic variables
syms x1 x2 x3 x4 x5 x6 x7 u1 u2 u3 u4 u5 C_av m R

% Define the matrix h(X)
H = [
    x1;
    x2;
    x3;
    (1/m)*((u2 + u3)*cos(u1)/R - (x4 + x5)*sin(u1) + (u4 + u5)/R - C_av*x1^2);
    (1/m)*((u2 + u3)*sin(u1)/R + (x4 + x5)*cos(u1) + x6 + x7)
];

% Define the state variables as a vector
state_vars = [x1; x2; x3; x4; x5; x6; x7];

% Compute the Jacobian matrix with respect to state variables
J_state = jacobian(H, state_vars);



% Display the Jacobian matrices
disp('The Jacobian matrix with respect to state variables is:');
disp(J_state);

