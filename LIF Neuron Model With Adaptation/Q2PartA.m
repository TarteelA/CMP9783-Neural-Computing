%Adaptive Exponential Leaky Integrate And Fire Model
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

%Voltage Range For Spike Uptick
Delta_th = 2.0e-3;

%Membrane Resistance (MΩ)
R_m = 100.0e6;

%Membrane Capacitance (uF/cm^2)
C_m = 100.0e-12;

%Potassium Reversal Potential
E_K = -80.0e-3;

%Leak Conductance (mS/cm^2)
gL = 10.0e-9;

%I_sra Control Term
a = 2.0e-9;

%I_sra Current Step
b = 2.0e-2;

%Time Constant For G_sra
T_sra = 0.2;

%Time Parameters
%Time Step (ms)
dt = 0.0001;

%Total Simulation Time (ms)
T = 1.5;

%Time Vector (ms)
t = 0:dt:T;  

%Applied Current Of 500 (pA)
I_app = zeros(1, length(t));
I_app(5001:10001) = 500e-12;

%Initialize Membrane Potential
V = zeros(1, length(t));
V(1) = E_L;     

%Spike Rate Adaptation Conductunce
I_sra = zeros(1, length(t));

%1.5s Neuron Simulation Loop
for i = 1:length(t)-1
    if (V(i) > V_th)
        V(i) = V_reset;                 
        I_sra(i) = I_sra(i) + b;
    end
        %Update Membrane Potential Using Euler Method
        V(i+1) = V(i) + dt * (gL * (E_L - V(i) + Delta_th * exp((V(i) - V_th)/Delta_th)) - I_sra(i) + I_app(i))/C_m;
        I_sra(i+1) = I_sra(i) + dt * (a * (V(i) - E_L) - I_sra(i))/T_sra;
end

figure(3);

%Plot Applied Current And Membrane Potential Over Time
subplot(2,1,1);
plot(t, 1e12*I_app);
xlabel('Time (s)');
ylabel('Applied Current (pA)');
title('Applied Current And Membrane Potential');

subplot(2,1,2);
plot(t, 1000*V);
xlabel('Time (s)');
ylabel('Membrane Potential (mV)');
