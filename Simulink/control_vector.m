function [avg_steer, Tfl, Tfr, Trl, Trr] = control_vector(left_steer, right_steer, Tdfl, Tdfr, Tdrl, Tdrr, Tbfl, Tbfr, Tbrl, Tbrr)
    
    driving_T = [Tdfl; Tdfr; Tdrl; Tdrr];
    braking_T = [Tbfl; Tbfr; Tbrl; Tbrr];
    total_T = driving_T - braking_T;
    
    avg_steer = (left_steer + right_steer) / 2;
    Tfl = total_T(1);
    Tfr = total_T(2);
    Trl = total_T(3);
    Trr = total_T(4);
    
end
