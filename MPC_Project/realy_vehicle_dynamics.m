function [Ux_dot, Uy_dot, r_dot] = real_vehicle_dynamics(Fxr, Fxf, Fzf, Fzr, delta, veh, state)
    muf = veh.muf;  mur = veh.mur;
    Caf = veh.Caf; Car = veh.Car;
    a = veh.a; b = veh.b;
    Ux = state.Ux; Uy = state.Uy; r = state.r;
    Iz = veh.Iz; m = veh.m;
    
    alpha_slidef = calc_sliding_angle(muf, Fzf, Caf);
    alpha_slider = calc_sliding_angle(mur, Fzr, Car);
    [alphaf, alphar] = calc_slip_angles(delta, Uy, a, b, r, Ux);

    Fyf = lateral_tire_force(alphaf, Caf, muf, mu_sf, Fzf, alpha_slidef);
    Fyr = lateral_tire_force(alphar, Car, mur, mu_sr, Fzr, alpha_slider);
    
    % equations of motion
    Ux_dot = (Fxr + Fxf*cos(delta) - Fyf*sin(delta) + m*r*Uy)/m;
    Uy_dot = (Fyf*cos(delta) + Fyr + Fxf*sin(delta) - m*r*Ux)/m;
    r_dot = (a*Fyf*cos(delta) + a*Fxf*sin(delta) - b*Fyr)/Iz;
end