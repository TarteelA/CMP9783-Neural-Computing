%Leaky Integrate And Fire Model With An Adaptation Current
%Part A: Simulate Mode Neuron For 1.5s With Current Pulse Of Iapp:500pA From 0.5s Until 1.0s
%Clear Previous Runs And Environment
clf;
clear;
close all;
clc;

%LIF Neuron Parameters
%Leak Reversal Potential (mV)
E_L = -75.0e-3;

%Threshold Potential (mV)
V_th = -50.0e-3;

%Reset Potential (mV)
V_reset = -80.0e-3;

%Membrane Resistance (MΩ)
R_m = 100.0e6;

%Membrane Capacitance (uF/cm^2)
C_m = 100.0e-12;

%Potassium Reversal Potential
E_K = -80.0e-3;

%Step Change In G_sra
Delta_G_sra = 1.0e-9;

%Time Constant For G_sra
T_sra = 0.2;

%Time Parameters
%Time Step (ms)
dt = 0.0001;

%Total Simulation Time (ms)
T = 1.5;

%Time Vector (ms)
t = 0:dt:T;  

%Initialize Membrane Potential
V = zeros(1, length(t));
V(1) = E_L;     

%Spike Rate Adaptation Conductunce
G_sra = zeros(1, length(t));

%Applied Current Of 500 (pA)
I_app = zeros(1, length(t));
I_app(5001:10001) = 500e-12;

%Simulation Loop
for i = 1:length(t)-1
   %Update Spike Handling In Simulation Loop:
    if (V(i) > V_th)
        %Set Spike Value For Visualization
        V(i) = V_reset;                 
        G_sra(i) = G_sra(i) + Delta_G_sra;
    end
        %Update Membrane Potential Using Euler Method
        V(i+1) = V(i) + dt * ((E_L - V(i))/R_m + G_sra(i) * (E_K - V(i)) + I_app(i))/C_m;
        G_sra(i+1) = G_sra(i) - dt * (G_sra(i)/T_sra);
end

figure(1);

%Plot Applied Current Over Time
subplot(3,1,1);
plot(t, 1e12*I_app);
xlabel('Time (ms)');
ylabel('I app (pA)');
title('Applied Current');

%Plot Membrane Potential Over Time
subplot(3,1,2);
plot(t, 1000*V);
xlabel('Time (ms)');
ylabel('V m (mV)');
title('Membrane Potential');

%Plot Adaptation Conductance Over Time
subplot(3,1,3);
plot(t, G_sra*1e9);
xlabel('Time (ms)');
ylabel('G sra (nS)');
title('Adaptation Conductance');