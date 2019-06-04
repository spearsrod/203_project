function alpha_sl = calc_sliding_angle(mu, Fz, Ca)
    alpha_sl = atan(3*mu*Fz/Ca);
end