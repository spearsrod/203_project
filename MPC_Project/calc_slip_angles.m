function [alphaf, alphar] = calc_slip_angles(delta, Uy, a, b, r, Ux)
    alphaf = (Uy + a*r)/Ux - delta;
    alphar = (Uy - b*r)/Ux;
end