%Author: Tarteel Alkaraan (25847208)
%Last Updated: 01/12/2024
%Leaky Integrate And Fire Model With Adaptation Current
%Part A: Simulate Mode Neuron For 1.5s With Current Pulse Of Iapp:500pA From 0.5s Until 1.0s
%Clear Previous Runs And Environment
clf;
clear;
close all;
clc;

%Time Parameters
%Time Step (ms)
dt = 0.0001;
%Total Simulation Time (ms)
T = 1.5;
%Time Vector (ms)
t = 0:dt:T;  

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
        V(i+1) = V(i) + dt * ((E_L - V(i)) / R_m + G_sra(i) * (E_K - V(i)) + I_app(i)) / C_m;
        G_sra(i+1) = G_sra(i) - dt * (G_sra(i) / T_sra);
end

figure(1);

%Plot Applied Current Over Time
subplot(3,1,1);
plot(t, 1e12 * I_app);
xlabel('Time (ms)');
ylabel('I app (pA)');
title('Applied Current');
grid on;

%Plot Membrane Potential Over Time
subplot(3,1,2);
plot(t, 1000 * V);
xlabel('Time (ms)');
ylabel('V m (mV)');
title('Membrane Potential');
grid on;

%Plot Adaptation Conductance Over Time
subplot(3,1,3);
plot(t, G_sra * 1e9);
xlabel('Time (ms)');
ylabel('G sra (nS)');
title('Adaptation Conductance');
grid on;
figure();

%Part B: Simulate Model For 5s With Range Of 20 Different Levels Of Constant Applied Current
%Time Parameters
%Time Step (ms)
dt = 0.0001;
%Total Simulation Time (ms)
T = 5;
%Time Vector (ms)
t = 0:dt:T; 

%Applied Current (pA)
I_app = 0:0.005:0.55;
I_app = I_app * 1e-9;

%Array For Storing 1/First ISI
Start_rate = zeros(size(I_app));

%Array For Storing 1/Last ISI
Last_rate = zeros(size(I_app));

%Array For Storing "One" For Just One Spike
Single_spike = zeros(size(I_app));

for j = 1:length(I_app)
    %Initialize Membrane Potential
    V = zeros(1, length(t));
    V(1) = E_L;  

    %Spike Rate Adaptation Conductunce
    G_sra = zeros(1, length(t));

    %Initialize Spikes
    Spikes = zeros(size(t));

    %Time Vector Simulation Loop
    for i = 1:length(t)-1
        if (V(i) > V_th)
            V(i) = V_reset;                 
            G_sra(i) = G_sra(i) + Delta_G_sra;
            Spikes(i) = 1;
        end
        %Update Membrane Potential Using Euler Method
        V(i+1) = V(i) + dt * ((E_L - V(i)) / R_m + G_sra(i) * (E_K - V(i)) + I_app(j)) / C_m;
        G_sra(i+1) = G_sra(i) - dt * (G_sra(i) / T_sra);
    end
    %Extract Spike Times
    Spike_times = dt * find(Spikes);

    if (length(Spike_times) > 1)
        %Interval Between Spikes
        ISI_s = diff(Spike_times);
    
        %Inverse Of First ISI
        Start_rate(j) = 1 / ISI_s(1);
        if (length(ISI_s) > 1)
            %Inverse Of Last ISI
            Last_rate(j) = 1 / ISI_s(end);
        end

    else
        if(isscalar(Spike_times))
            Single_spike(j) = 1;
        end
    end
end

hold on;

%Plot Applied Current And Spike Rate
plot(I_app * 1e9, Last_rate);

ISI_indes = find(Start_rate);
plot(1e9 * I_app(ISI_indes), Start_rate(ISI_indes), 'o');

ISI_indes = find(Single_spike);
plot(1e9 * I_app(ISI_indes), 0 * Single_spike(ISI_indes), 'sk');

xlabel('I app (nA)');
ylabel('Spike Rate (Hz)');

legend('Last Rate', '1/ISI(1)', 'Single Spike');

title('Applied Current And Spike Rate');
grid on;

hold off;
figure(2);