function path = monza_path()

factor = 70;
s_path = [0 5 6.5 8 10 25 29 31 33 37 42 45 48 59 62 65 68 75 81 87];
k_path = [0 -1/1.5 1 0 -1/9 0 1/2 -1/4 0 -1/5 0 -1/7.5 0 1/2.5 -1/3 1/5 0 -1/5 0 0];
s_path = s_path*factor;
k_path = k_path/factor;
% integrate s and k with given initial conditions [psi0, E0, N0] to get path
path = generate_path(s_path,k_path,[0;0;0]);
end