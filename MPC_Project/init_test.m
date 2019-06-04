clear all; close all;
veh.Caf = 188000; veh.Car = 203000;
veh.m = 1600; veh.Iz = 2200;
veh.L = 2.4;
veh.wf_percent = 0.577; veh.wr_percent = 0.423;
g = 9.81;
a = veh.L*veh.wr_percent; b = veh.L*veh.wf_percent;
veh.a = a;
veh.b = b;
veh.rW = 0.35;
veh.K_phi_f = 90000;
veh.K_phi_r = 65000;
veh.COG = [0 500];

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

N = 10;
K = 1/20;
K_dot = 0;
Ux = 10;
h = 0.01;
A = [0 0 0 0 0; 0 0 1 0 0; 0 0 -(Caf + Car)/(m*Ux) (Caf + Car)/m (b*Car - a*Caf)/(m*Ux); 0 0 0 0 1; 0 0 (b*Car - a*Caf)/(Iz*Ux) (a*Caf - b*Car)/Iz -(a^2*Caf + b^2*Car)/(Iz*Ux)];
A = A*h + eye(5);
B = [0; 0; Caf/m; 0; a*Caf/Iz];
B = B*h;

C = [Ux; 0; -K*(Ux^2 - (b*Car-a*Caf)/m); 0; -K*(a^2*Caf + b^2*Car)/Iz - K_dot*Ux];
C = C*h;
C = repelem(C, 1, N);


umin = deg2rad(-25); umax = deg2rad(25);
e_max = 2; e_min = -2;
e_dot_max = e_max*10; e_dot_min = e_min*10;
dpsi_max = pi/3; dpsi_min = -pi/3;
dpsi_dot_max = 10*dpsi_max; dpsi_dot_min = 10*dpsi_min;
xmax = [100; e_max; e_dot_max; dpsi_max; dpsi_dot_max];
xmin = [-101; e_min; e_dot_min; dpsi_min; dpsi_dot_min];

path_length = 100;
x = [-path_length; 0; 0; 0; 0];
%cvx precision
cvx_precision(min(max(min(abs(x))/(10),1e-6),0.9999))

cvx_begin quiet
    variables X(5,N+1) U(1,N)
    max(X') <= xmax'; max(U') <= umax';
    min(X') >= xmin'; min(U') >= umin';
    size(C)
    X(:,2:N+1) == A*X(:,1:N)+B*U + C;
    %X(:,2:N+1) == A*X(:,1:N)+B*U;
    X(:,1) == x; %initial state constraint
%     if(final == 1)
%         X(:,N+1) == 0;
%     end
%     cost1 = 0;
%     cost2 = 0;
%     for j = 1:N    X    X
%         cost1 = cost1 + X(:,j).'*Q*X(:,j);
%         cost2 = cost2 + U(j)*R*U(j);
%     end
%    J = X(:,N+1).'*P*X(:,N+1) + cost1 + cost2;
    J = X(1,N+1).'*X(1,N+1);
    minimize(J)
cvx_end

cvx_status
if strcmp(cvx_status,'Solved')

    %store controlUx = 10;
    u = U

    %accumulate cost
    optvalmpc = X(1,N+1)

    %forward propagate state
    x = X

    %record state
    Xallmpc = x

else
   % break from loop
   optvalmpc = Inf
end