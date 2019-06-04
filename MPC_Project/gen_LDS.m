function [A, B, C] = gen_LDS(K, K_dot, veh, Ux, h)
a = veh.a; b = veh.b;
Caf = veh.Caf; Car = veh.Car;
Iz = veh.Iz; m = veh.m;
A = [0 K*Ux 0 Ux 0; 0 0 1 0 0; 0 0 -(Caf + Car)/(m*Ux) (Caf + Car)/m (b*Car - a*Caf)/(m*Ux); 0 0 0 0 1; 0 0 (b*Car - a*Caf)/(Iz*Ux) (a*Caf - b*Car)/Iz -(a^2*Caf + b^2*Car)/(Iz*Ux)];
A = A*h + eye(5);
B = [0; 0; Caf/m; 0; a*Caf/Iz];
B = B*h;

C = [Ux; 0; -K*(Ux^2 - (b*Car-a*Caf)/m); 0; -K*(a^2*Caf + b^2*Car)/Iz - K_dot*Ux];
C = C*h;
%C = repelem(C, 1, N);
end