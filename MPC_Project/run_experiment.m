function [e_m, e_h] = run_experiment(veh, path_center, t_final, dT, Ux_des, e_0)
%Set vehicle specs
Caf = veh.Caf; Car = veh.Car;
m = veh.m; Iz = veh.Iz;
L = veh.L;
wf_percent = veh.wf_percent; wr_percent = veh.wr_percent;
g = 9.81;
wf = m*g*wf_percent; wr = m*g*wr_percent;
a = L*wr_percent; b = L*wf_percent;
muf = 0.97; mu_sf = 0.97;
mur = 1.03; mu_sr = 1.03;
Fzf = wf; Fzr = wr;
a = veh.a;
b = veh.b;
rW = veh.rW;

t_s = 0:dT:t_final;
N = length(t_s);

% allocate space for simulation data
dpsi_rad    = zeros(N,1);
s_m         = zeros(N,1);
e_m         = zeros(N,1);
Ux_m        = zeros(N,1);
Uy_m        = zeros(N,1);
r_m         = zeros(N,1);
delta_rad   = zeros(N,1);
e_h = zeros(N,1);
%e_lh = zeros(N,1);


% set initial conditions
dpsi_rad(1) = 0;
s_m(1)      = 0;
e_m(1)      = e_0;
Ux_m(1)     = Ux_des;
Uy_m(1)     = 0;
r_m(1)      = 0;



K_grad = (wf/Caf - wr/Car)/g;
% simulation loop
for idx = 1:N
    % look up K
    %s_m(idx)
    K = interp1(path.s_m, path.k_1pm, s_m(idx));

    % current states
    dpsi = dpsi_rad(idx);
    s = s_m(idx);
    e = e_m(idx);
    Ux = Ux_m(idx);
    Uy = Uy_m(idx);
    r = r_m(idx);
    
    % compute delta
    Kla = 3500;
    xla = 15;
    Cafd = 188000;
    Carff = 203000;
    dpsi_ss = K*((m*a*Ux^2)/(L*Carff) - b);
    delta_ff = Kla*xla/Cafd*dpsi_ss + K*(L + K_grad*Ux^2);
    delta = -Kla*(e + xla*dpsi)/Cafd + delta_ff;
%    delta = -Kla*(e + xla*dpsi)/Cafd;
    %TODO
    %Change this so it can bump on a single wheel (average delta change)
    if(bump_profile(idx,1) ~= 0)
        bump_size = bump_profile(idx,1);
        bump_side = bump_profile(idx,2);
        if(bump_side == 0)
            [dtheta, final_right_suspension] = BumpSteer(init_suspension.right_sus, bump_size);
            final_suspension.left_sus = init_suspension.left_sus;
            final_suspension.right_sus = final_right_suspension;
        else
            [dtheta, final_left_suspension] = BumpSteer(init_suspension.left_sus, bump_size);
            final_suspension.right_sus = init_suspension.right_sus;
            final_suspension.left_sus = final_left_suspension;
        end
        delta = delta + dtheta/2;
        if(~isreal(dtheta))
            e_m = 10000*ones(N,1);
            e_h = 10000*ones(N,1);
            return
        end
    end
    
    
    %Compute Fx's
    Kd = m*0.1*9.81;
    Fxtotal = Kd*(Ux_des - Ux);
    Fxf = 0.5*Fxtotal; Fxr = Fxtotal*0.5;
    if(K == 0)
        alpha_slidef = calc_sliding_angle(muf, Fzf, Caf);
        alpha_slider = calc_sliding_angle(mur, Fzr, Car);
        [alphaf, alphar] = calc_slip_angles(delta, Uy, a, b, r, Ux);
        
        Fyf = lateral_tire_force(alphaf, Caf, muf, mu_sf, Fzf, alpha_slidef);
        Fyr = lateral_tire_force(alphar, Car, mur, mu_sr, Fzr, alpha_slider);
    else
        velocity = sqrt(Ux^2 + Uy^2);
        [Fzif, Fzof, Fzir, Fzor] = roll_displacement(init_suspension, veh, velocity, K, Fzf, Fzr);
        % TODO, make sure the tire stiffness is reasonable when splitting like
        % this
        alpha_slideif = calc_sliding_angle(muf, Fzif, Caf);
        alpha_slideir = calc_sliding_angle(mur, Fzir, Car);
        alpha_slideof = calc_sliding_angle(muf, Fzof, Caf);
        alpha_slideor = calc_sliding_angle(mur, Fzor, Car);
        [alphaf, alphar] = calc_slip_angles(delta, Uy, a, b, r, Ux);


        Fyif = lateral_tire_force(alphaf, Caf, muf, mu_sf, Fzif, alpha_slideif);
        Fyir = lateral_tire_force(alphar, Car, mur, mu_sr, Fzir, alpha_slideir);
        Fyof = lateral_tire_force(alphaf, Caf, muf, mu_sf, Fzof, alpha_slideof);
        Fyor = lateral_tire_force(alphar, Car, mur, mu_sr, Fzor, alpha_slideor);
        
        Fyf = (Fyif + Fyof)/2;
        Fyr = (Fyir + Fyor)/2;
        
    end
    
    % equations of motion
    Ux_dot = (Fxr + Fxf*cos(delta) - Fyf*sin(delta) + m*r*Uy)/m;
    Uy_dot = (Fyf*cos(delta) + Fyr + Fxf*sin(delta) - m*r*Ux)/m;
    r_dot = (a*Fyf*cos(delta) + a*Fxf*sin(delta) - b*Fyr)/Iz;
    
    %r_dot
    
    s_dot = (1/(1 - e*K))*(Ux*cos(dpsi) - Uy*sin(dpsi));
    e_dot = Uy*cos(dpsi) + Ux*sin(dpsi);
    dpsi_dot = r - K*s_dot;


    % only update next state if we are not awf/2t end of simulation
    delta_rad(idx) = delta;
    if idx < N
        % euler integration
        dpsi_rad(idx+1) = dpsi_rad(idx) + dpsi_dot*dT;
        s_m(idx+1)      = s_m(idx) + s_dot*dT;
        e_m(idx+1)      = e_m(idx) + e_dot*dT; 
        Ux_m(idx+1)     = Ux_m(idx) + Ux_dot*dT;
        Uy_m(idx+1)     = Uy_m(idx) + Uy_dot*dT;
        r_m(idx+1)      = r_m(idx) + r_dot*dT;
    end
    
    e_h(idx) = xla*dpsi; 
end
%animate(path, veh, dpsi_rad, s_m, e_m, delta_rad)

end