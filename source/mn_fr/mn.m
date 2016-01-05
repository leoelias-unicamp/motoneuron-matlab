function v = mn(t,V)
%MN Differential equations and dynamic of motoneuron model.
% V = MN(T,V) returns the somatic and dendritic membrane potentials after
% integration.

global SPIKE0 t0 SPIKE1 t1 pd rp %Monitor index_monitor
global gc El Cd gld Cs gls Vth 
global ECa GCa ENa GNa EK GKf GKs
global m h n q p m0 h0 n0 q0 p0 
global alfa_M beta_M alfa_H beta_H alfa_N beta_N alfa_Q beta_Q alfa_P beta_P
global re we Ge Esyn_exc
% global ri wi Gi Esyn_ini

%% Pulse-based constant rates
if V(2) > Vth && (t0 == -1 || (t-t0) > rp)
    SPIKE0 = 1;
    t0 = t;
    m0 = m;
    h0 = h;
    n0 = n;
    q0 = q;
end
if SPIKE0 == 1 && (t-t0) > pd
    SPIKE0 = 0;
    t0 = t;
    m0 = m;
    h0 = h;
    n0 = n;
    q0 = q;
end

if SPIKE0 > 0
    m = 1+(m0-1)*exp(-alfa_M*(t-t0));
    h = h0*exp(-beta_H*(t-t0));
    n = 1+(n0-1)*exp(-alfa_N*(t-t0));
    q = 1+(q0-1)*exp(-alfa_Q*(t-t0));
else
    m = m0*exp(-beta_M*(t-t0));
    h = 1+(h0-1)*exp(-alfa_H*(t-t0));
    n = n0*exp(-beta_N*(t-t0));
    q = q0*exp(-beta_Q*(t-t0));
end

% if V(2) > Vth-3.00 && (t1 == -1 || (t-t1) > 0.05)
%     SPIKE1 = 1;
%     t1 = t;
%     p0 = p;
% end
% if SPIKE1 == 1 && (t-t1) > pd
%     SPIKE1 = 0;
%     t1 = t;
%     p0 = p;
% end
% 
% if SPIKE1 > 0
%     p = 1+(p0-1)*exp(-alfa_P*(t-t1)); 
% else
%     p = p0*exp(-beta_P*(t-t1)); 
% end

%% Synaptic current
% Excitatory
re(1) = r_exc(1,we(1,:),t);
re(2) = r_exc(2,we(2,:),t);
re(3) = r_exc(3,we(3,:),t);
re(4) = r_exc(4,we(4,:),t);
% re(5) = r_exc(5,we(5,:),t);
% re(6) = r_exc(6,we(6,:),t);
% re(7) = r_exc(7,we(7,:),t);
% re(8) = r_exc(8,we(8,:),t);
% re(9) = r_exc(9,we(9,:),t);
% re(10) = r_exc(10,we(10,:),t);
% re(11) = r_exc(11,we(11,:),t);
% re(12) = r_exc(12,we(12,:),t);
% re(13) = r_exc(13,we(13,:),t);
% re(14) = r_exc(14,we(14,:),t);
% re(15) = r_exc(15,we(15,:),t);
% re(16) = r_exc(16,we(16,:),t);
% re(17) = r_exc(17,we(17,:),t);
% re(18) = r_exc(18,we(18,:),t);
% re(19) = r_exc(19,we(19,:),t);
% re(20) = r_exc(20,we(20,:),t);
% re(21) = r_exc(21,we(21,:),t);
% re(22) = r_exc(22,we(22,:),t);
% re(23) = r_exc(23,we(23,:),t);
% re(24) = r_exc(24,we(24,:),t);
% re(25) = r_exc(25,we(25,:),t);
% re(26) = r_exc(26,we(26,:),t);
% re(27) = r_exc(27,we(27,:),t);
% re(28) = r_exc(28,we(28,:),t);
% re(29) = r_exc(29,we(29,:),t);
% re(30) = r_exc(30,we(30,:),t);
% re(31) = r_exc(31,we(31,:),t);
% re(32) = r_exc(32,we(32,:),t);
% re(33) = r_exc(33,we(33,:),t);
% re(34) = r_exc(34,we(34,:),t);
% re(35) = r_exc(35,we(35,:),t);
% re(36) = r_exc(36,we(36,:),t);
% re(37) = r_exc(37,we(37,:),t);
% re(38) = r_exc(38,we(38,:),t);
% re(39) = r_exc(39,we(39,:),t);
% re(40) = r_exc(40,we(40,:),t);
% re(41) = r_exc(41,we(41,:),t);
% re(42) = r_exc(42,we(42,:),t);
% re(43) = r_exc(43,we(43,:),t);
% re(44) = r_exc(44,we(44,:),t);
% re(45) = r_exc(45,we(45,:),t);
% re(46) = r_exc(46,we(46,:),t);
% re(47) = r_exc(47,we(47,:),t);
% re(48) = r_exc(48,we(48,:),t);
% re(49) = r_exc(49,we(49,:),t);
% re(50) = r_exc(50,we(50,:),t);
% re(51) = r_exc(51,we(51,:),t);
% re(52) = r_exc(52,we(52,:),t);
% re(53) = r_exc(53,we(53,:),t);
% re(54) = r_exc(54,we(54,:),t);

%Inhibitory
% ri(1) = r_ini(1,wi(1,:),t);
% ri(2) = r_ini(2,wi(2,:),t);
% ri(3) = r_ini(3,wi(3,:),t);
% ri(4) = r_ini(4,wi(4,:),t);
% ri(5) = r_ini(5,wi(5,:),t);
% ri(6) = r_ini(6,wi(6,:),t);
% ri(7) = r_ini(7,wi(7,:),t);
% ri(8) = r_ini(8,wi(8,:),t);
% ri(9) = r_ini(9,wi(9,:),t);
% ri(10) = r_ini(10,wi(10,:),t);
% ri(11) = r_ini(11,wi(11,:),t);
% ri(12) = r_ini(12,wi(12,:),t);
% ri(13) = r_ini(13,wi(13,:),t);
% ri(14) = r_ini(14,wi(14,:),t);
% ri(15) = r_ini(15,wi(15,:),t);
% ri(16) = r_ini(16,wi(16,:),t);
% ri(17) = r_ini(17,wi(17,:),t);
% ri(18) = r_ini(18,wi(18,:),t);
% ri(19) = r_ini(19,wi(19,:),t);
% ri(20) = r_ini(20,wi(20,:),t);
% ri(21) = r_ini(21,wi(21,:),t);
% ri(22) = r_ini(22,wi(22,:),t);
% ri(23) = r_ini(23,wi(23,:),t);
% ri(24) = r_ini(24,wi(24,:),t);
% ri(25) = r_ini(25,wi(25,:),t);
% ri(26) = r_ini(26,wi(26,:),t);
% ri(27) = r_ini(27,wi(27,:),t);

%% Noise
% noise = -0.9+var_syn.*randn(1,1);

%% Integration
% Dendritic Membrane Potential - default
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)+Iinj_d(t));

% Dendritic Membrane Potential - passive dendrite
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)+Iinj_d(t));

% Dendritic Membrane Potential - constant Ge
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)-Ge*(V(1)-Esyn_exc)+Iinj_d(t));
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-Ge*(V(1)-Esyn_exc)+Iinj_d(t));

% Dendritic Membrane Potential - current noise
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)+Iinj_d(t)+noise);

% Dendritic Membrane Potential - random activity of synaptic inputs, only excitatory
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)-Ge*sum(re)*(V(1)-Esyn_exc)+Iinj_d(t));
v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-Ge*sum(re)*(V(1)-Esyn_exc)+Iinj_d(t));

% Dendritic Membrane Potential - random activity of synaptic inputs, only inhibitory
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)-Gi*sum(ri)*(V(1)-Esyn_ini)+Iinj_d(t));

% Dendritic Membrane Potential - random activity of synaptic inputs, excitatory and inhibitory inputs
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)-Ge*sum(re)*(V(1)-Esyn_exc)-Gi*sum(ri)*(V(1)-Esyn_ini)+Iinj_d(t));

% Somatic Membrane Potential
% In voltage-clamp protocols uncomment (B), and comment (A). If you want to
% simulate voltage ramp set v(2,1) as the slope value, and Vs0 (main.m) as
% 0.

% (A):
v(2,1) = (1/Cs)*(-gc*(V(2)-V(1))-gls*(V(2)-El)-GNa*m^3*h*(V(2)-ENa)-GKf*n^4*(V(2)-EK)-GKs*q^2*(V(2)-EK)+Iinj_s(t));

% (B):
% Ramp
% if t < 6000
%     v(2,1) = 0.008;                                                      
% else
%     v(2,1) = -0.008;
% end

% Constant
% v(2,1) = 0.0;

%% Monitoring Conductances
% Monitor(index_monitor,1) = GNa*m^3*h;
% Monitor(index_monitor,2) = GKf*n^4;
% Monitor(index_monitor,3) = GKs*q^2;
% Monitor(index_monitor,4) = GCa*p;
% index_monitor = index_monitor+1;