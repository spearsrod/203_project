


A = [0 1 0 0; 0 -(Caf + Car)/(m*Ux) (Caf + Car)/m (b*Car - a*Caf)/(m*Ux); 0 0 0 1; 0 (b*Car - a*Caf)/(Iz*Ux) (a*Caf - b*Car)/Iz -(a^2*Caf + b^2*Car)/(Iz*Ux)];
B = [0; Caf/m; 0; a*Caf/Ix];
C = [0; -K*(Ux^2 - (b*Car-a*Caf)/m); 0; -K*(a^2*Caf + B^2*Car)/Iz - K_dot*Ux];


%Goal minimize -s(T)

