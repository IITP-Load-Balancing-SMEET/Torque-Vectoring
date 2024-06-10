% Define the symbolic variables
syms x1 x2 x3 x4 x5 x6 x7 u1 u6 u7 u8 u9 C_av m

% Define the matrix h(X)
H = [
    x1;
    x2;
    x3;
    (1/m)*((u6 + u7)*cos(u1) - (x4 + x5)*sin(u1) + u8 + u9 - C_av*x1^2);
    (1/m)*((u6 + u7)*sin(u1) + (x4 + x5)*cos(u1) + x6 + x7)
];

% Define the state variables as a vector
state_vars = [x1; x2; x3; x4; x5; x6; x7];

% Compute the Jacobian matrix with respect to state variables
J_state = jacobian(H, state_vars);



% Display the Jacobian matrices
disp('The Jacobian matrix with respect to state variables is:');
disp(J_state);

