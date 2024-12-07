%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Adaptive Leaky Integrate-and-Fire (LIF) Neuron Model
% This script simulates the LIF neuron with adaptation and analyzes its behavior.
% Part A: Simulates the effect of a 500 pA pulse on the neuron for 1.5s.
% Part B: Examines firing behavior across a range of applied currents and computes the f-I curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
clear all; % Clear all variables from the workspace
close all; % Close all open figures

% Parameters for LIF model
E_l = -75e-3;         % Leak reversal potential (V)
V_th = -50e-3;        % Threshold potential (V)
V_reset = -80e-3;     % Reset potential (V)
R_m = 100e6;          % Membrane resistance (Ω)
C_m = 100e-12;        % Membrane capacitance (F)
E_k = -80e-3;         % Adaptation reversal potential (V)
Delta_GSRA = 1e-9;    % Spike-triggered adaptation increment (S)
Tau_SRA = 0.2;        % Adaptation time constant (s)

%% Simulation parameters for Part A
dt = 0.0001;          % Time step (s)
T = 1.5;              % Total simulation time (s)
time = 0:dt:T;        % Time vector

% Input current: 500 pA pulse from 0.5s to 1.0s
I_app = zeros(size(time)); % Initialize applied current as zero
I_app(time >= 0.5 & time <= 1.0) = 500e-12; % Apply 500 pA during the specified interval

% Initialize variables
V = E_l * ones(size(time));  % Membrane potential initialized to leak potential
G_sra = zeros(size(time));   % Adaptation conductance initialized to 0

% Simulation loop for Part A
for i = 1:length(time)-1
    % Check if the membrane potential exceeds the threshold
    if V(i) > V_th
        % Reset membrane potential and increase adaptation conductance if a
        % spike is generated
        V(i) = V_reset;                
        G_sra(i) = G_sra(i) + Delta_GSRA; % Increment of adaptation conductance
    end
    
    % Update membrane potential and adaptation conductance
    %Calculate currents
    I_leak = (E_l - V(i)) / R_m;       % Leak current due to passive ion flow
    I_adapt = G_sra(i) * (E_k - V(i)); % Adaptation current opposing spiking activity

    % Updating membrane potential using Euler's method
    dVdt = (I_leak + I_adapt + I_app(i)) / C_m; % Total current's effect on membrane potential
    V(i+1) = V(i) + dVdt * dt; % Euler integration step for membrane potential

    % Update adaptation conductance with decay
    G_sra(i+1) = G_sra(i) - dt * (G_sra(i) / Tau_SRA); % Decay term
    G_sra(i+1) = max(0, G_sra(i+1));  % Ensure conductance stays non-negative
end

% Plot results for Part A
figure;
subplot(3, 1, 1);
plot(time, I_app * 1e12);
ylabel('I_{app} (pA)', 'FontWeight', 'normal');
title('Input Current vs. Time', 'FontWeight', 'normal');
grid on;

subplot(3, 1, 2);
plot(time, V * 1e3);
ylabel('V_m (mV)', 'FontWeight', 'normal');
title('Membrane Potential vs. Time', 'FontWeight', 'normal');
grid on;

subplot(3, 1, 3);
plot(time, G_sra * 1e9);
xlabel('Time (s)', 'FontWeight', 'normal');
ylabel('G_{sra} (nS)', 'FontWeight', 'normal');
title('Adaptation Conductance vs. Time', 'FontWeight', 'normal');
grid on;

sgtitle('LIF Model with Adaptation', 'FontSize', 16, 'FontWeight', 'bold');

%% Section 1, Question 1B:

%% Simulation parameters for Part B
T_B = 5;                % Total simulation time (s)
time_B = 0:dt:T_B;      % Time vector
I_range = linspace(240e-12, 550e-12, 20);  % Range of applied currents

% Initialize storage for results
initial_rate = zeros(size(I_range)); % Initial firing rate
steady_rate = zeros(size(I_range));  % Steady-state firing rate
single_spike = zeros(size(I_range)); % Tracks instances where only one spike occurs

% Simulation loop for Part B
for j = 1:length(I_range)
    % Reset variables for each input current level
    V = E_l * ones(size(time_B));    % Reset membrane potential
    G_sra = zeros(size(time_B));    % Reset adaptation conductance
    spikes = zeros(size(time_B));   % Initialize spike tracker as 0

    for i = 1:length(time_B)-1
        % Check for spike
        if V(i) > V_th
             % Reset potential and increment adaptation conductance if a
             % spike is detected
            V(i) = V_reset;          % Reset potential
            G_sra(i) = G_sra(i) + Delta_GSRA; % Increment adaptation conductance
            spikes(i) = 1;          % Record spike
        end

        % Calculate currents
        I_leak = (E_l - V(i)) / R_m; % Leak current
        I_adapt = G_sra(i) * (E_k - V(i)); % Adaptation current
        dVdt = (I_leak + I_adapt + I_range(j)) / C_m; % Total current's effect
        V(i+1) = V(i) + dVdt * dt; % Update membrane potential
        G_sra(i+1) = G_sra(i) - dt * (G_sra(i) / Tau_SRA); % Decay adaptation conductance
        G_sra(i+1) = max(0, G_sra(i+1));  % Ensure conductance stays non-negative
    end

    % Analyze spike times
    spike_times = find(spikes) * dt; % Spike times in seconds
    if length(spike_times) > 1
        ISIs = diff(spike_times);         % Interspike intervals
        initial_rate(j) = 1 / ISIs(1);    % First ISI
        steady_rate(j) = 1 / mean(ISIs(end-5:end)); % Compute steady-state rate
    elseif length(spike_times) == 1
        single_spike(j) = 1;              % Record single spike occurrence
    end
end

% Plot f-I Curve for Part B
figure;
hold on;
plot(I_range * 1e12, steady_rate, 'b-', 'LineWidth', 2, 'DisplayName', 'Steady-State Firing Rate');
plot(I_range * 1e12, initial_rate, 'ro', 'MarkerSize', 6, 'LineWidth', 1.5, 'DisplayName', '1 / Initial ISI');
plot(I_range(single_spike == 1) * 1e12, zeros(sum(single_spike), 1), 'g*', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Single Spike');
xlabel('I_{app} (pA)', 'FontWeight', 'bold');
ylabel('Firing Rate (Hz)', 'FontWeight', 'bold');
title('LIF f-I Curve (Part B)', 'FontWeight', 'bold');
legend('show');
grid on;
hold off;