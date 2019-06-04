clear all; close all;


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

