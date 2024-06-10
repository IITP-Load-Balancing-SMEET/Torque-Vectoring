function [X, P] = update(z, Xp, Pp, u)
    persistent Cav

    % H is 5 X 11 matrix
    H = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
         0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0;
         -2*Cav/m, 0, 0, cos(u(1))/m, cos(u(1))/m, 1/m, 1/m, -sin(u(1))/m, -sin(u(1))/m, 0, 0;
         0, 0, 0, sin(u(1))/m, sin(u(1))/m, 0, 0, cos(u(1))/m, cos(u(1))/m, 1/m, 1/m];

    K = (Pp * H') / (H * Pp * H' + R);
    err = z - h(Xp, u);
    X = Xp + (K * err);
    P = Pp - (K * H * Pp);
end

function z_hat = h(Xp, u)
    persistent Cav

    z_hat = zeros(5, 1);
    z_hat(1:4) = Xp(1:4);

    z_hat(4) = 1/m * (((Xp(4) + Xp(5)) * cos(u(1))) - ...
        ((Xp(8) + Xp(9)) * sin(u(1))) + ...
        (Xp(6)) + ...
        (Xp(7)) - ...
        (Cav * x(1)^2));

    z_hat(5) = 1/m * (((Xp(4) + Xp(5)) * sin(u(1))) - ...
        ((Xp(8) + Xp(9)) * cos(u(1))) + ...
        (Xp(10)) + ...
        (Xp(11)));

end