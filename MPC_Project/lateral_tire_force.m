function Fy = lateral_tire_force(alpha, Ca, mu, mu_s, Fz, alpha_slide)
    if(abs(alpha) < alpha_slide)
        Fy = -Ca*tan(alpha) + Ca^2/(3*mu*Fz)*(2 - mu_s/mu)*abs(tan(alpha))*tan(alpha) - Ca^3/(9*mu^2*Fz^2)*tan(alpha)^3*(1 - 2*mu_s/(3*mu));
    else
        Fy = -mu_s*Fz*sign(alpha);
    end
end