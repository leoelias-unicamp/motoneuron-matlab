%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% S-type motoneuron model with active dendrite %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
clear all
clc

tic
global s SPIKE0 SPIKE1 t0 t1 pd rp Monitor index_monitor
global gc El Cd gld Cs gls Vth 
global ECa GCa ENa GNa EK GKf GKs
global m h n q p m0 h0 n0 q0 p0
global alfa_M beta_M alfa_H beta_H alfa_N beta_N alfa_Q beta_Q alfa_P beta_P
% global re time_inst1 Ge Esyn_exc
% global te0 de fe j0 X_exc SYNAPSE_exc t2 pd1 r0_exc aux0
% global r_inf_exc tau_R_exc beta_R_exc
% global ri time_inst2 Gi Esyn_ini
% global ti0 di j1 X_ini SYNAPSE_ini t3 pd2 r0_ini aux1
% global r_inf_ini tau_R_ini beta_R_ini

% global var_freq

% load W.mat                                                                 %Load ISI from a file

%% Simulation Time (ms)
ts = 0.0;                                                                  %Start Time
s = 0.01;                                                                  %Step
te = 10000.0;                                                               %End Time
t = ts:s:te;

%% Monitoring Conductances
Monitor = zeros(4*length(t),4);
index_monitor = 1;

%% Model Parameters
Irh = 5.0;                                                                 %Rheobase Current [nA]
Ri = 70.0;                                                                 %Cytoplasm Resistivity [Ohm x cm]
Cm = 1.0;                                                                  %Membrane Specific Capacitance [uF/cm?]
Rm_d = 12500.0;                                                            %Dendritic Membrane Specific Resistance [Ohm x cm?]
Rm_s = 1100.0;                                                             %Somatic Membrane Specific Resistance [Ohm x cm?]
ld = 6154.0;                                                               %Dendritic Length [um]
ls = 80;                                                                   %Somatic Length [um]
rd = (52.0)/2;                                                             %Dendritic Radius [um]
rs = (80.0)/2;                                                             %Somatic Radius [um]
pd = 0.6;                                                                  %Pulse Duration [ms]
rp = 5.0;                                                                  %Refractory Period [ms]
SPIKE0 = 0.0;
t0 = -1.0;
t1 = -1.0;
SPIKE1 = 0.0;

%% Common Constants
gc = (2*1e2)/(((Ri*ld)/(pi*(rd^2)))+((Ri*ls)/(pi*(rs^2))));                %Coupling Conductance [uS]
El = 0.0;                                                                  %Leakage Voltage [mV]

%% Dendritic-compartment Constants
Cd = 2*1e-5*pi*rd*ld*Cm;                                                   %Dendritic Capacitance [nF]
gld = (2*1e-2*pi*rd*ld)/Rm_d;                                              %Leakage Conductance [uS]
p = 0;                                                                     %Activation variable for L-type Calcium Conductance
ECa = 140;                                                                 %Equilibrium Potential - L-type Calcium [mV]
GCa = 2*1e-5*pi*rd*ld*(0.038);                                             %Peak Conductance - L-type Calcium [uS]

%% Somatic-compartment Constants
Cs = 2*1e-5*pi*rs*ls*Cm;                                                   %Somatic Capacitance [nF]
gls = (2*1e-2*pi*rs*ls)/Rm_s;                                              %Leakage Conductance [uS]
m = 0;                                                                     %Activation Variable for Sodium Conductance
h = 1;                                                                     %Inactivation Variable for Sodium Conductance
ENa = 120;                                                                 %Equilibrium Potential - Sodium [mV]
GNa = 2*1e-5*pi*rs*ls*(30.0);                                              %Peak Conductance - Sodium [uS] //MAX: 30.0
EK = -10;                                                                  %Equilibrium Potential - Potassium [mV]
n = 0;                                                                     %Activation Variable for Fast Potassium Conductance
q = 0;                                                                     %Activation Variable for Slow Potassium Conductance
GKf = 2*1e-5*pi*rs*ls*(4.0);                                               %Peak Conductance - Fast Potassium [uS] //MAX: 4.0
GKs = 2*1e-5*pi*rs*ls*(19.15);                                             %Peak Conductance - Slow Potassium [uS] //MAX: 19.15
Rn = 1/(gls+((gld*gc)/(gld+gc)));                                          %Input Resistance [Mega Ohm]
Vth = floor(Rn*Irh);                                                       %Membrane Potential Threshold [mV]

%% Ionic Channels Constants
alfa_M = 22.00;                                                            %[1/ms]
beta_M = 13.00;                                                            %[1/ms]
alfa_H = 0.500;                                                            %[1/ms]
beta_H = 4.000;                                                            %[1/ms]
alfa_N = 1.500;                                                            %[1/ms]
beta_N = 0.100;                                                            %[1/ms]
alfa_Q = 1.500;                                                            %[1/ms]
beta_Q = 0.025;                                                            %[1/ms]
alfa_P = 0.008;                                                            %[1/ms]
beta_P = 0.014;                                                            %[1/ms]

%% Synaptic Current Constants
% Tmax = 1.0;                                                                %Transmitter [mM]

% Excitatory
% te0 = 500.0;
% AFe = 1;                                                                   %Afferent Ia fibers //MAX: 54
% Ge = 600.0e-3;                                                             %Peak Conductance - Excitatory Synapse [uS]
% Esyn_exc = 70;                                                             %Reversal Potential - Excitatory Synapse [mV]
% alfa_R_exc = 0.5;                                                          %[1/ms]
% beta_R_exc = 2.5;                                                          %[1/ms]
% r_inf_exc = (alfa_R_exc*Tmax)/(alfa_R_exc*Tmax+beta_R_exc);
% tau_R_exc = 1/(alfa_R_exc*Tmax+beta_R_exc);
% pd1 = 0.2;                                                                 %Pulse Duration [ms]
% re = zeros(AFe,1);
% X_exc = zeros(AFe,1);
% SYNAPSE_exc = zeros(AFe,1);
% t2 = -1.0*ones(AFe,1);
% aux0 = zeros(AFe,1);
% j0 = ones(AFe,1);

% Inhibitory
% ti0 = 0.0;                  
% AFi = 1;                                                                   %Afferent IN-Ia fibers //MAX: 27
% Gi = 500.0e-3;                                                             %Peak Conductance - Inhibitory Synapse [uS]
% Esyn_ini = -16;                                                            %Reversal Potential - Inhibitory Synapse [mV]
% alfa_R_ini = 0.5;                                                          %[1/ms]
% beta_R_ini = 1.2;                                                          %[1/ms]
% r_inf_ini = (alfa_R_ini*Tmax)/(alfa_R_ini*Tmax+beta_R_ini);
% tau_R_ini = 1/(alfa_R_ini*Tmax+beta_R_ini);
% pd2 = 0.4;                                                                 %Pulse Duration [ms]
% ri = zeros(AFi,1);
% X_ini = zeros(AFi,1);
% SYNAPSE_ini = zeros(AFi,1);
% t3 = -1.0*ones(AFi,1);
% aux1 = zeros(AFi,1);
% j1 = ones(AFi,1);

%% ISI [ms]
% Excitatory
% lambda_exc = 0.0;                                                        %Pre-synaptic neurons firing rate [pps]

% Poisson
% for i = 1:AFe
%     for j = 1:length(t)
%         we(i,j) = -(1000/lambda_exc)*log(1-rand(1,1));
%         if we(i,j) < 0.2
%             we(i,j) = 0.2;
%         end
%         if j == 1
%             time_inst1(i,j) = te0 + we(i,j);
%         else
%             time_inst1(i,j) = time_inst1(i,j-1) + we(i,j);
%         end
%         if time_inst1(i,j) > te
%             break
%         end
%     end
% end

% Sinusoidal Modulation
% Fmod = var_freq;                                                           %Frequency of the modulation [Hz]
% TimeAux = 0;
% for i = 1:AFe
%     de(i+1,1) = (i)*(1000/(lambda_exc*AFe));
%     for j = 1:length(t)
%         fe(1,j) = lambda_exc + 0.20.*lambda_exc.*sin(2.*pi.*(Fmod/1000).*TimeAux);
%         we(1,j) = 1000./fe(1,j);
%         if we(1,j) < 0.2
%             we(1,j) = 0.2;
%         end
%         if j == 1
%             time_inst1(i,j) = de(i,j);
%             TimeAux = time_inst1(1,j) + we(1,j);
%         else
%             de(i+1,j) = (i)*(we(1,j)/AFe);
%             time_inst1(1,j) = time_inst1(1,j-1) + we(1,j);
%             TimeAux = time_inst1(1,j);
%             time_inst1(i+1,j) = time_inst1(1,j) + de(i+1,j);
%         end
%         if time_inst1(i,j) > te
%             break
%         end
%     end
% end
% 
% time_inst1 = time_inst1(1:AFe,:);

% Inhibitory
% lambda_ini = 0.0;                                                          %Pre-synaptic neurons firin rate [pps]

% Poisson
% for i = 1:AFe
%     for j = 1:length(t)
%         wi(i,j) = -(1000/lambda_ini)*log(1-rand(1,1));
%         if wi(i,j) < 0.2
%             wi(i,j) = 0.2;
%         end
%         if j == 1
%             time_inst2(i,j) = te0 + wi(i,j);
%         else
%             time_inst2(i,j) = time_inst2(i,j-1) + wi(i,j);
%         end
%         if time_inst2(i,j) > te
%             break
%         end
%     end
% end

%% Numerical Integration - 4th Order Runge-Kuta
% Initial Conditions
m0 = 0.0;                                                                  %Activation Variable for Sodium Conductance
h0 = 1.0;                                                                  %Inactivation Variable for Sodium Conductance
n0 = 0.0;                                                                  %Activation Variable for Fast Potassium Conductance
q0 = 0.0;                                                                  %Activation Variable for Slow Potassium Conductance
p0 = 0.0;                                                                  %Activation Variable for L-type Calcium Conductance
% r0_exc = zeros(AFe,1);                                                     %Fraction of Bound Receptors of Excitatory Synapse
% r0_ini = zeros(AFi,1);                                                     %Fraction of Bound Receptors of Inhibitory Synapse
Vd0 = 0.0;                                                                 %Dendrite Membrane Potential
Vs0 = 0.0;                                                                 %Somatic Membrane Potential
IC = [Vd0; Vs0];

% In voltage-clamp protocols set Vs0 as a constant value. Comment (A), and uncomment (B) in the function 'mn.m' 

V = ode4(@mn,t,IC);

%% Resampling Ionic Conductances
Monitor = downsample(Monitor,4);

%% Figures
figure
plot(t,V(:,1),'k')
grid
title('Dendritic Membrane Potential')
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')

figure
plot(t,V(:,2),'k')
grid
title('Somatic Membrane Potential')
xlabel('Time (ms)')
ylabel('Membrane Potential (mV)')

% In = gc.*(V(:,2)-V(:,1));
% figure
% plot(t,In,'k')
% grid
% xlabel('Membrane Potential (mV)')
% ylabel('Membrane Current (nA)')
% Ileak = 0.4427.*V(:,2)+0.0227;
% Ieff = In-Ileak;
% Ipic = max(-Ieff(1:600000))
% Isus = max(-Ieff(600001:end))

%% Firing Rate
finst

%% Twitch
% spike = zeros(1,length(t));
% Apeak = 1.70;                                                              %Peak twitch [gf]
% tpeak = 89.0;                                                              %Contraction time [ms] 
% Tetanic = 10.00; 
% 
% for z=1:length(index)
%     spike(1,index) = 1/s;
% end
% 
% % Filter coefficients 
% a(1) = 1;
% a(2) = -2*exp(-(s/tpeak));
% a(3) = exp(-2*(s/tpeak));
% b(1) = 0;
% b(2) = Apeak*((s^2)/tpeak)*exp(1-(s/tpeak));
% 
% Tw = filter(b,a,spike);
% 
% for zz=1:length(Tw)
%     if Tw(zz) > Tetanic
%         Tw(zz) = Tetanic;
%     end
% end
% 
% figure
% plot(t,Tw,'k')
% grid
% title('Motor Unit Twitch')
% xlabel('Time (ms)')
% ylabel('Force (gf)')

toc