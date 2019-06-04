function [s_dot, e_dot, dpsi_dot] = path_dynamics(path_state, veh_state, K)
    %r_dot
    
    s = path_state.s; e = path_state.e; dpsi = path_state.dpsi;
    Ux = veh_state.Ux; Uy = veh_state.Uy; r = veh_state.r;
    
    s_dot = (1/(1 - e*K))*(Ux*cos(dpsi) - Uy*sin(dpsi));
    e_dot = Uy*cos(dpsi) + Ux*sin(dpsi);
    dpsi_dot = r - K*s_dot;
end