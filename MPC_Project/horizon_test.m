clear all; close all;
%% Simple Finite Horizon Test
veh.Caf = 188000; veh.Car = 203000;
veh.m = 1600; veh.Iz = 2200;
veh.L = 2.4;
veh.wf_percent = 0.577; veh.wr_percent = 0.423;
g = 9.81;
a = veh.L*veh.wr_percent; b = veh.L*veh.wf_percent;
veh.a = a;
veh.b = b;
veh.rW = 0.35;

Fxmax = 4*veh.m;
Fxmin = 4*veh.m;


Caf = veh.Caf; Car = veh.Car;
m = veh.m; Iz = veh.Iz;
L = veh.L;
wf_percent = veh.wf_percent; wr_percent = veh.wr_percent;
wf = m*g*wf_percent; wr = m*g*wr_percent;
a = L*wr_percent; b = L*wf_percent;
muf = 0.97; mu_sf = 0.97;
mur = 1.03; mu_sr = 1.03;
veh.muf = muf; veh.mur = mur;
Fzf = wf; Fzr = wr;
a = veh.a;
b = veh.b;
rW = veh.rW;

T_s = 300;
N = 100;
K = 1/20;
%Ux = 10;
K_dot = 0;

h = 0.05;
% A = [0 0 0 0 0; 0 0 1 0 0; 0 0 -(Caf + Car)/(m*Ux) (Caf + Car)/m (b*Car - a*Caf)/(m*Ux); 0 0 0 0 1; 0 0 (b*Car - a*Caf)/(Iz*Ux) (a*Caf - b*Car)/Iz -(a^2*Caf + b^2*Car)/(Iz*Ux)];
% A = A*h + eye(5);
% B = [0; 0; Caf/m; 0; a*Caf/Iz];
% B = B*h;
% 
% C = [Ux; 0; -K*(Ux^2 - (b*Car-a*Caf)/m); 0; -K*(a^2*Caf + b^2*Car)/Iz - K_dot*Ux];
% C = C*h;
% C = repelem(C, 1, N);

path_length = 300;
path_length
umin = deg2rad(-25); umax = deg2rad(25);
e_max = 2; e_min = -2;
e_dot_max = e_max*10; e_dot_min = e_min*10;
dpsi_max = pi/3; dpsi_min = -pi/3;
dpsi_dot_max = 10*dpsi_max; dpsi_dot_min = 10*dpsi_min;
xmax = [path_length; e_max; e_dot_max; dpsi_max; dpsi_dot_max];
xmin = [-path_length; e_min; e_dot_min; dpsi_min; dpsi_dot_min];


s_path = [0 300 500 700];
k_path = [-K 0 K 0];

% integrate s and k with given initial conditions [psi0, E0, N0] to get path
path = generate_path(s_path,k_path,[0;0;0]);
%path = monza_path();


s0 = -path_length; e0 = 0; dpsi0 = 0;
x0 = [s0; e0; 0; dpsi0; 0];
path_state.s = s0 + path_length;
path_state.e = e0; path_state.dpsi = dpsi0;
x = x0;

Ux0 = 0.1; Uy0 = 0; r0 = 0;
veh_state.Ux = Ux0; veh_state.Uy = Uy0; veh_state.r = r0;


optvalmpc = 0;

%store solutionspath_length = h*Ux*T_s;
Xallmpc = zeros(5,T_s+1); Uallmpc = zeros(1,T_s);
Xallmpc(:,1) = x;
s_profile = zeros(N,1);
for idx = 2:N
    s_profile(idx) = s_profile(idx - 1) + 10*h;
end
%K_prev = path.k_1pm(1);
for idx = 1:T_s
    %cvx precision
    idx
    K_prev = interp1(path.s_m, path.k_1pm, path_state.s);
    cvx_precision(min(max(min(abs(x))/(10),1e-6),0.9999))

    Q = 1/(idx);
    
    cvx_begin quiet
        variables Uxdes(1,N+1) Fx(1,N) Sdes(1,N)
        max(Fx') <= Fxmax;
        min(Fx') >= Fxmin; 
        cost = 0;
        for j = 1:N
            s_cur = s_profile(j);
            K_cur = interp1(path.s_m, path.k_1pm, s_cur);
            R = abs(1/(K_cur + 0.00001))
            Uxmax = sqrt(R*muf*(Fzf + Fzr))
            max(Uxdes(j)) <= Uxmax;
            min(Uxdes(j)) >= 0;
            Uxdes(j+1) == Uxdes(j)+1/m*Fx(j)*h;
            cost = cost + Uxdes(j+1);
        end
        Uxdes(1) == veh_state.Ux; %initial state constraint
        J = -cost;
        minimize(J);
    cvx_end

    if strcmp(cvx_status,'Solved')
        s_des = zeros(N,1);
        s_des(1) = path_state.s;
        for idx = 2:N
            s_des(idx) = s_des(idx - 1) + Uxdes(idx)*h;
        end
    else
       % break from loop
       hi = 2
       optvalmpc = Inf;
       break;
    end
    Uxdes
    s_des
    cvx_begin quiet
        variables X(5,N+1) U(1,N)
        max(X') <= xmax'; max(U') <= umax';
        min(X') >= xmin'; min(U') >= umin';
        cost2 = 0;
        for j = 1:N
            s_cur = s_des(j);
            Ux = Uxdes(j);
            K_cur = interp1(path.s_m, path.k_1pm, s_cur);
            K_dot_cur = K_cur - K_prev;
            [A, B, C] = gen_LDS(K_cur, K_dot_cur, veh, Ux, h);
            X(:,j+1) == A*X(:,j)+B*U(j)+C;
            K_prev = K_cur;
            %cost2 = cost2 + X(2,j+1).'*Q*X(2,j+1);
        end
        %X(:,2:N+1) == A*X(:,1:N)+B*U + C;
        %X(:,2:N+1) == A*X(:,1:N)+B*U;
        X(:,1) == x; %initial state constraint
        J = X(1,N+1).'*X(1,N+1);% + cost2;
        minimize(J);
    cvx_end
    cvx_status
    if strcmp(cvx_status,'Solved')

        %store control
        delta = U(1);
        Fxr = 0;
        Fxf = 0;
        
        [Ux_dot, Uy_dot, r_dot] = real_vehicle_dynamics(Fxr, Fxf, Fzf, Fzr, delta, veh, veh_state);
        [s_dot, e_dot, dpsi_dot] = path_dynamics(path_state, veh_state, K);
        
        Ux_new = veh_state.Ux + Ux_dot*h;
        Uy_new = veh_state.Uy + Uy_dot*h;
        r_new = veh_state.r + r_dot*h;
        s_new = path_state.s + s_dot*h;
        e_new = path_state.e + e_dot*h;
        dpsi_new = path_state.dpsi + dpsi_dot*h;
        
        veh_state.Ux = Ux_new; veh_state.Uy = Uy_new; veh_state.r = r_new;
        path_state.s = s_new; path_state.e = e_new; path_state.dpsi = dpsi_new;
        optimalvalmpc = s_new;
        x_model = A*x + B*delta + C;
        x = [s_new - path_length; e_new; x_model(3); dpsi_new; x_model(5)];
        Xallmpc(:,idx) = x;
        Uallmpc(idx) = delta;
        s_profile = X(1,2:N+1) + path_length;
    else
       % break from loop
       hi = 0
       optvalmpc = Inf;
       break;
    end
end
s_m = Xallmpc(1,:) + path_length;
e_m = Xallmpc(2,:);
dpsi_rad = Xallmpc(4,:);
delta_rad = Uallmpc(:);
animate(path, veh, dpsi_rad, s_m, e_m, delta_rad)