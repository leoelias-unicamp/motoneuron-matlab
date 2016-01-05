function v = mn(t,V)
%MN Differential equations and dynamic of motoneuron model.
% V = MN(T,V) returns the somatic and dendritic membrane potentials after
% integration.

global SPIKE0 SPIKE1 t0 t1 pd rp Monitor index_monitor
global gc El Cd gld Cs gls Vth 
global ECa GCa ENa GNa EK GKf GKs
global m h n q p m0 h0 n0 q0 p0 
global alfa_M beta_M alfa_H beta_H alfa_N beta_N alfa_Q beta_Q alfa_P beta_P
% global re time_inst1 Ge Esyn_exc
% global ri time_inst2 Gi Esyn_ini

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

% if V(2) > Vth-5.25 && (t1 == -1 || (t-t1) > 0.05)
% if V(2) > 2.75 && (t1 == -1 || (t-t1) > 0.05)
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
% re(1) = r_exc(1,time_inst1(1,:),t);
% re(2) = r_exc(2,time_inst1(2,:),t);
% re(3) = r_exc(3,time_inst1(3,:),t);
% re(4) = r_exc(4,time_inst1(4,:),t);
% re(5) = r_exc(5,time_inst1(5,:),t);
% re(6) = r_exc(6,time_inst1(6,:),t);
% re(7) = r_exc(7,time_inst1(7,:),t);
% re(8) = r_exc(8,time_inst1(8,:),t);
% re(9) = r_exc(9,time_inst1(9,:),t);
% re(10) = r_exc(10,time_inst1(10,:),t);
% re(11) = r_exc(11,time_inst1(11,:),t);
% re(12) = r_exc(12,time_inst1(12,:),t);
% re(13) = r_exc(13,time_inst1(13,:),t);
% re(14) = r_exc(14,time_inst1(14,:),t);
% re(15) = r_exc(15,time_inst1(15,:),t);
% re(16) = r_exc(16,time_inst1(16,:),t);
% re(17) = r_exc(17,time_inst1(17,:),t);
% re(18) = r_exc(18,time_inst1(18,:),t);
% re(19) = r_exc(19,time_inst1(19,:),t);
% re(20) = r_exc(20,time_inst1(20,:),t);
% re(21) = r_exc(21,time_inst1(21,:),t);
% re(22) = r_exc(22,time_inst1(22,:),t);
% re(23) = r_exc(23,time_inst1(23,:),t);
% re(24) = r_exc(24,time_inst1(24,:),t);
% re(25) = r_exc(25,time_inst1(25,:),t);
% re(26) = r_exc(26,time_inst1(26,:),t);
% re(27) = r_exc(27,time_inst1(27,:),t);
% re(28) = r_exc(28,time_inst1(28,:),t);
% re(29) = r_exc(29,time_inst1(29,:),t);
% re(30) = r_exc(30,time_inst1(30,:),t);
% re(31) = r_exc(31,time_inst1(31,:),t);
% re(32) = r_exc(32,time_inst1(32,:),t);
% re(33) = r_exc(33,time_inst1(33,:),t);
% re(34) = r_exc(34,time_inst1(34,:),t);
% re(35) = r_exc(35,time_inst1(35,:),t);
% re(36) = r_exc(36,time_inst1(36,:),t);
% re(37) = r_exc(37,time_inst1(37,:),t);
% re(38) = r_exc(38,time_inst1(38,:),t);
% re(39) = r_exc(39,time_inst1(39,:),t);
% re(40) = r_exc(40,time_inst1(40,:),t);
% re(41) = r_exc(41,time_inst1(41,:),t);
% re(42) = r_exc(42,time_inst1(42,:),t);
% re(43) = r_exc(43,time_inst1(43,:),t);
% re(44) = r_exc(44,time_inst1(44,:),t);
% re(45) = r_exc(45,time_inst1(45,:),t);
% re(46) = r_exc(46,time_inst1(46,:),t);
% re(47) = r_exc(47,time_inst1(47,:),t);
% re(48) = r_exc(48,time_inst1(48,:),t);
% re(49) = r_exc(49,time_inst1(49,:),t);
% re(50) = r_exc(50,time_inst1(50,:),t);
% re(51) = r_exc(51,time_inst1(51,:),t);
% re(52) = r_exc(52,time_inst1(52,:),t);
% re(53) = r_exc(53,time_inst1(53,:),t);
% re(54) = r_exc(54,time_inst1(54,:),t);

%Inhibitory
% ri(1) = r_ini(1,time_inst2(1,:),t);
% ri(2) = r_ini(2,time_inst2(2,:),t);
% ri(3) = r_ini(3,time_inst2(3,:),t);
% ri(4) = r_ini(4,time_inst2(4,:),t);
% ri(5) = r_ini(5,time_inst2(5,:),t);
% ri(6) = r_ini(6,time_inst2(6,:),t);
% ri(7) = r_ini(7,time_inst2(7,:),t);
% ri(8) = r_ini(8,time_inst2(8,:),t);
% ri(9) = r_ini(9,time_inst2(9,:),t);
% ri(10) = r_ini(10,time_inst2(10,:),t);
% ri(11) = r_ini(11,time_inst2(11,:),t);
% ri(12) = r_ini(12,time_inst2(12,:),t);
% ri(13) = r_ini(13,time_inst2(13,:),t);
% ri(14) = r_ini(14,time_inst2(14,:),t);
% ri(15) = r_ini(15,time_inst2(15,:),t);
% ri(16) = r_ini(16,time_inst2(16,:),t);
% ri(17) = r_ini(17,time_inst2(17,:),t);
% ri(18) = r_ini(18,time_inst2(18,:),t);
% ri(19) = r_ini(19,time_inst2(19,:),t);
% ri(20) = r_ini(20,time_inst2(20,:),t);
% ri(21) = r_ini(21,time_inst2(21,:),t);
% ri(22) = r_ini(22,time_inst2(22,:),t);
% ri(23) = r_ini(23,time_inst2(23,:),t);
% ri(24) = r_ini(24,time_inst2(24,:),t);
% ri(25) = r_ini(25,time_inst2(25,:),t);
% ri(26) = r_ini(26,time_inst2(26,:),t);
% ri(27) = r_ini(27,time_inst2(27,:),t);

%% Noise
% noise = -0.9+60.*randn(1,1);

%% Integration
% Dendritic Membrane Potential - default
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-0.6*GCa*p*(V(1)-ECa)+Iinj_d(t));

% Dendritic Membrane Potential - passive dendrite
v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)+Iinj_d(t));

% Dendritic Membrane Potential - constant Ge
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)-Ge*(V(1)-Esyn_exc)+Iinj_d(t));

% Dendritic Membrane Potential - current noise
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)+Iinj_d(t)+noise);

% Dendritic Membrane Potential - random activity of synaptic inputs, only excitatory
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-GCa*p*(V(1)-ECa)-Ge*sum(re)*(V(1)-Esyn_exc)+Iinj_d(t));
% v(1,1) = (1/Cd)*(-gc*(V(1)-V(2))-gld*(V(1)-El)-Ge*sum(re)*(V(1)-Esyn_exc)+Iinj_d(t)); %Passive dendrite

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
Monitor(index_monitor,1) = GNa*m^3*h;
Monitor(index_monitor,2) = GKf*n^4;
Monitor(index_monitor,3) = GKs*q^2;
% Monitor(index_monitor,4) = GCa*p;
index_monitor = index_monitor+1;