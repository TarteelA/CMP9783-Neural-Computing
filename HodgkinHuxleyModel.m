%Hodgkin Huxley Model As An Oscillator
%Clear Previous Runs And Environment
clf;
clear;
close all;
clc;

%Choices Are Question Tasks 2,3,4,5, And 6
Task = '2';

%Time Parameters
%Time Step (ms)
dt = 2e-8;  
%Total Simulation Time (ms)
T = 0.35;    
%Time Vector
time = 0:dt:T;  

%Hodgkin-Huxley Parameters
%Membrane Capacitance (uF/cm^2)
C_m = 100e-12;   

%Maximum Conductance Of Sodium (mS/cm^2)
g_Na = 12e-6; 

%Maximum Conductance Of Potassium (mS/cm^2)
g_K = 3.6e-6;   

%Leak Conductance (mS/cm^2)
g_L = 30e-9;    

%Sodium Reversal Potential (mV)
E_Na = 0.045;    

%Potassium Reversal Potential (mV)
E_K = -0.082;    

%Leak Reversal Potential (mV)
V_L = -0.060;  

%Initial Time Applied Current Begins
I_begin = 100e-3; 

%Initial Length Of Applied Current Pulse
I_Len = 5e-3; 

%Initial Baseline Current Before/After Pulse
I_Base = 0e-9; 

%Initial Number Of Current Pulses
Num_Pulses = 1; 

%Initial Separation Between Current Pulses
Sep_Pulse = 20e-3; 

%Default Value For V
V_zero = -0.065;

%Default Value For M
M_zero = 0.05;

%Default Value For H
H_zero = 0.5;

%Default Value For N
N_zero = 0.35;

%Edit Few Parameters Or Default Conditions With Task Part
switch Task
    case '2'
        %Length Of Applied Current Pulse In Task 2
        I_Len = 100e-3;
        %Applied Current During Pulse In Task 2
        I_E = 0.22e-9;
    case '3'
        %10 Pulses In Task 3
        Num_Pulses = 10;
        %Applied Current For Every Pulse In Task 3 (Like As 2)
        I_E = 0.22e-9;
        %Time Between Pulses In Task 3
        Sep_Pulse = 18e-3;
    case '4'
        %10 Pulses In Task 4
        Num_Pulses = 10;
        %Large Baseline Current In Task 4
        I_Base = 0.6e-9;
        %Pulse Current Is Below Baseline In Task 4
        I_E = 0e-9;
    case '5'
        %Baseline Current In Task 5
        I_Base = 0.65e-9;
        %Pulsed Current In Task 5
        I_E = 1e-9;
    case '6'
        %Baseline Current In Task 6
        I_Base = 0.7e-9;
        %Pulsed Current In Part 6 (Like As 5)
        I_E = 1e-9;
        %Default Condition For M In Task 6
        M_zero = 0;
        %Default Condition For H In Task 6
        H_zero = 0;
        %Default Condition For N In Task 6
        N_zero = 0;
end

%Declare Applied Current Vector
%Define Current Vector At Baseline
I_app = I_Base * ones(size(time));

%For Every Current Pulse
for Pulse = 1:Num_Pulses
    %Time Onset Of Pulse
    Pulse_Begin = I_begin + (Pulse - 1) * Sep_Pulse;
    %Time Offset Of Pulse
    Pulse_End = Pulse_Begin + I_Len;
    %Define Applied Current Value For Duration Of Current Pulse
    for I = round(Pulse_Begin / dt) + 1:round(Pulse_End / dt)
        I_app(I) = I_E;
    end
end

%Initialize Variables
%Membrane Potential (mV)
V = zeros(size(time));
V(1) = V_zero;

%Sodium Activation
M = zeros(size(time));
M(1) = M_zero;

%Sodium Inactivation
H = zeros(size(time));
H(1) = H_zero;

%Potassium Activation
N = zeros(size(time));
N(1) = N_zero;

%Total Current
I_Tot = zeros(size(time));

%Sodium Current
I_Na = zeros(size(time));

%Potassium Current
I_K = zeros(size(time));
I_L = zeros(size(time));

%Simulation Loop
for I = 2:length(time)

    V_M = V(I-1);
    
    %Define Rate Functions For Gating Variables
    if (V_M == -0.045)
        A_M = 1e3;
    else
        A_M = (1e5 * (- V_M - 0.045)) / (exp(100 * (- V_M - 0.045)) -1);
    end
    B_M = 4000 * exp((- V_M - 0.070) / 0.018);
    A_H = 70 * exp(50 * (- V_M - 0.070));
    B_H = 1000 / (1 + exp(100 * (- V_M - 0.040)));

    if (V_M == -0.060)
        A_N = 100;
    else
        A_N = (1e4 * (- V_M - 0.060)) / (exp(100 * (- V_M - 0.060)) -1);
    end
    B_N = 125 * exp((- V_M - 0.070) / 0.08);
    
    %Sodium Activation Gating Variable
    tau_M = 1 / (A_M + B_M);
    M_INF = A_M / (A_M + B_M);   
    
    %Sodium Inactivation Gating Variable
    tau_H = 1 / (A_H + B_H);
    H_INF = A_H / (A_H + B_H);    
    
    %Potassium Activation Gating Variable
    tau_N = 1 / (A_N + B_N);
    N_INF = A_N / (A_N + B_N);   

    %Update Gating Variables Using Euler Method
    M(I) = M(I-1) + (M_INF - M(I-1)) * dt / tau_M;
    H(I) = H(I-1) + (H_INF - H(I-1)) * dt / tau_H;
    N(I) = N(I-1) + (N_INF - N(I-1)) * dt / tau_N;
    
    %Compute Ionic Currents
    I_Na(I) = g_Na * M(I) * M(I) * M(I) * H(I) * (E_Na - V(I-1));
    I_K(I) = g_K * N(I) * N(I) * N(I) * N(I) * (E_K - V(I-1));
    I_L(I) = g_L * (V_L - V(I-1));
    I_Tot(I) = I_L(I) + I_Na(I) + I_K(I) + I_app(I);
    
    %Update Membrane Potential Using Euler's Method
    V(I) = V(I-1) + I_Tot(I) * dt / C_m;
end

%Plot Applied Current And Membrane Potential Over Time
hold on;

figure(1)

subplot(2,1,1)
plot(time(10:10:end), 1e9 * I_app(10:10:end))
Title_Str = ['Task ', Task, newline, ' Applied Current'];
title(Title_Str)
xlabel('Time (ms)')
ylabel('I app (nA)')
axis([0 T -0.5 1.05])
hold on

subplot(2,1,2)
plot(time(10:10:end), 1e3 * V(10:10:end))
Title_Str2 = 'Membrane Potential';
title(Title_Str2)
xlabel('Time (ms)')
ylabel('V m (mV)')

if (max(V) > 0)
    axis([0 T -85 45])
else
    axis([0 T -80 -55])
end