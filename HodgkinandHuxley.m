%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hodgkin-Huxley Model Simulation
% This script simulates the Hodgkin-Huxley (HH) neuron model with various 
% applied current stimuli and analyzes the resulting membrane potential and gating variables.
% It includes multiple parts (A-F), each demonstrating different input scenarios.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Define global parameters
global G_leak G_na G_k E_na E_k V_leak C_m dt
G_leak = 30e-9;
G_na = 12e-6;
G_k = 3.6e-6;
E_na = 45e-3;
E_k = -82e-3;
V_leak = -60e-3;
C_m = 100e-12;
dt = 0.0001e-3; % Time step in seconds

% Setup time vector
t = 0:dt:0.35; % Total simulation time

%% Part A: Stabilization at -70.2 mV
% Simulates stabilization of the HH model at -70.2 mV with no applied current.

% Applied current
I_app = zeros(1, length(t));

% Initial Conditions
v_init = -70.2e-3; % Initial membrane potential (-70.2 mV)
m_init = 0;        % Initial sodium activation variable
h_init = 0;        % Initial sodium inactivation variable
n_init = 0;        % Initial potassium activation variable

% Run the HH simulation
[v_sim, m_sim, h_sim, n_sim] = hhsim(t, I_app, v_init, m_init, h_init, n_init);

% Generate plots for Part A
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app * 1e9); % Convert current to nA for readability
xlabel("Time (s)", 'FontSize', 12); 
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part A: Applied Current", 'FontSize', 14, 'FontWeight', 'bold');
grid on; 

% Plot membrane potential
subplot(2,1,2);
plot(t, v_sim * 1e3); % Convert voltage to mV for readability
xlabel("Time (s)", 'FontSize', 12); 
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part A: Membrane Potential", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

%% Part B: Subthreshold Oscillation
% Simulates the HH model with a 0.22 nA pulse for 100 ms starting at 100 ms.

% Applied current for Part B
I_app((t >= 0.1) & (t < 0.2)) = 0.22e-9; % 0.22 nA pulse for 100 ms starting at 100 ms

% Run simulation with initial conditions for Part B
v_init = -61e-3; % Default initial condition

% Run the HH simulation
[v_sim, m_sim, h_sim, n_sim] = hhsim(t, I_app, v_init, m_init, h_init, n_init);

% Generate plots for Part B
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app * 1e9); % Convert current to nA for readability
xlabel("Time (s)", 'FontSize', 12);
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part B: Applied Current", 'FontSize', 14, 'FontWeight', 'bold');
grid on; 

% Plot membrane potential
subplot(2,1,2);
plot(t, v_sim * 1e3); % Convert voltage to mV for readability
xlabel("Time (s)", 'FontSize', 12);
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part B: Membrane Potential", 'FontSize', 14, 'FontWeight', 'bold');
grid on; 

%% Part C: Repeated Pulses with 16 ms and 18 ms Delays
% Simulates the HH model with repeated pulses spaced 16 ms or 18 ms apart.

% Parameters for 16 ms delay
delay_16 = 16e-3;      % Delay between pulses (s)
pulse_duration = 5e-3; % Duration of each pulse (s)

% Create applied current vector for 16 ms delay.
I_app_16 = zeros(1, length(t));
delay_index_16 = floor(delay_16/dt); % Index corresponding to 16 ms
pulse_width = floor(pulse_duration/dt); % Index width for a 5 ms pulse

% Generate current pulses for 16 ms delay
for i = 1:10
    left = (i-1)*delay_index_16 + 1; % Start index of pulse
    right = min(left + pulse_width - 1, length(t)); % End index of pulse
    I_app_16(left:right) = 0.22e-9; % Apply current
end

% Run simulation for 16 ms delay.
[v_sim_16, m_sim_16, h_sim_16, n_sim_16] = hhsim(t, I_app_16);

% Plot results for 16 ms delay.
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app_16 * 1e9); % Convert current to nA
xlabel("Time (s)", 'FontSize', 12);
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part C: Applied Current vs. Time (Delay = 16 ms)", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot Membrane Potential 
subplot(2,1,2);
plot(t, v_sim_16 * 1e3); % Convert voltage to mV
xlabel("Time (s)", 'FontSize', 12);
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part C: Membrane Potential vs. Time (Delay = 16 ms)", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Parameters for 18 ms delay
delay_18 = 18e-3; % Delay between pulses (s)

% Create applied current vector for 18 ms delay
I_app_18 = zeros(1, length(t));
delay_index_18 = floor(delay_18/dt); % Index corresponding to 18 ms

% Generate current pulses for 18 ms delay
for i = 1:10
    left = (i-1)*delay_index_18 + 1; % Start index of pulse
    right = min(left + pulse_width - 1, length(t)); % End index of pulse
    I_app_18(left:right) = 0.22e-9; % Apply current
end

% Run simulation for 18 ms delay.
[v_sim_18, m_sim_18, h_sim_18, n_sim_18] = hhsim(t, I_app_18);

% Plot results for 18 ms delay.
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app_18 * 1e9); % Convert current to nA
xlabel("Time (s)", 'FontSize', 12);
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part C: Applied Current vs. Time (Delay = 18 ms)", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot Membrane Potential
subplot(2,1,2);
plot(t, v_sim_18 * 1e3); % Convert voltage to mV
xlabel("Time (s)", 'FontSize', 12);
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part C: Membrane Potential vs. Time (Delay = 18 ms)", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

%% Part D: Baseline Current with Periodic Dropouts
% Simulates the HH model with baseline current of 0.6 nA and periodic current dropouts.

% Create vector for applied current with baseline 0.6 nA.
I_app = zeros(1, length(t)) + 0.6e-9; % Baseline current (0.6 nA)
delay = 20e-3; % 20 ms delay between dropouts (converted to econds)
pulse_duration = 5e-3; % 5 ms pulse duration (converted to seconds)

% Calculate delay and pulse indices for current dropouts
delay_index = floor(delay / dt); % No of steps for delay
pulse_width = floor(pulse_duration / dt); % No of steps for pulse duration

% Set applied current to 0 during dropout intervals.
for i = 1:10
    left = (i-1)*delay_index + 1; % Start index for current dropout
    right = min(left + pulse_width - 1, length(t)); % End index for current dropout
    I_app(left:right) = 0; % Drop current to 0 during the pulse interval
end

% Run simulation with proper initial conditions.
v_init = -65e-3; % Initial membrane potential in volts (-65 mV)
m_init = 0.05;   % Initial value for 'm' (Sodium activation)
h_init = 0.6;    % Initial value for 'h' (Sodium inactivation)
n_init = 0.35;   % Initial value for 'n' (Potassium activation)

[v_sim, m_sim, h_sim, n_sim] = hhsim(t, I_app, v_init, m_init, h_init, n_init);

% Generate plots for part d.
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app * 1e9); % Convert current to nA for plotting
xlabel("Time (s)", 'FontSize', 12);
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part D: Applied Current with Periodic Dropouts", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot Membrane Potential
subplot(2,1,2);
plot(t, v_sim * 1e3); % Convert membrane potential to mV for plotting
xlabel("Time (s)", 'FontSize', 12);
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part D: Membrane Potential with Periodic Dropouts", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

%% Part E: Single Stimulus with Higher Amplitude
% Simulates the HH model with baseline current of 0.65 nA and a 1 nA pulse at 100 ms.

% Create vector for applied current with baseline 0.65 nA
I_app = zeros(1, length(t)) + 0.65e-9; % Baseline current

% Define stimulus parameters.
pulse_start_time = 100e-3; % Pulse starts at 100 ms (converted to seconds)
pulse_duration = 5e-3; % Pulse lasts 5 ms (converted to seconds)
pulse_amplitude = 1e-9; % Pulse amplitude 

% Calculate indices for the pulse.
pulse_start_index = floor(pulse_start_time / dt); % Start index of the pulse
pulse_end_index = pulse_start_index + floor(pulse_duration / dt); % End index of the pulse

% Apply the pulse.
I_app(pulse_start_index:pulse_end_index) = pulse_amplitude;

% Run simulation with realistic initial conditions
v_init = -65e-3; % Initial membrane potential in volts (-65 mV).
m_init = 0.05;   % Initial value for 'm' (typically small).
h_init = 0.6;    % Initial value for 'h' (typically around 0.6-0.8).
n_init = 0.35;   % Initial value for 'n' (typically around 0.3-0.5).

[v_sim, m_sim, h_sim, n_sim] = hhsim(t, I_app, v_init, m_init, h_init, n_init);

% Generate plots for applied current and membrane potential.
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app * 1e9); % Convert current to nA for plotting.
xlabel("Time (s)", 'FontSize', 12);
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part E: Applied Current with Single Stimulus", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot Membrane Potential
subplot(2,1,2);
plot(t, v_sim * 1e3); % Convert membrane potential to mV for plotting.
xlabel("Time (s)", 'FontSize', 12);
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part E: Membrane Potential with Single Stimulus", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Generate plots for gating variables.
figure;
plot(t, m_sim, 'LineWidth', 1.5);
hold on;
plot(t, h_sim, '-', 'LineWidth', 1.5);
plot(t, n_sim, '--', 'LineWidth', 1.5);
xlabel("Time(s)", 'FontSize', 12);
ylabel("Gating Variables", 'FontSize', 12);
title("Gating Variables vs. Time (Part E)", 'FontSize', 14, 'FontWeight', 'bold');
legend("Sodium Activation (m)", "Sodium Inactivation (h)", "Potassium Activation (n)");
grid on;

%% Part F: Single Stimulus with Different Initial Conditions
% Simulates the HH model with baseline current of 0.7 nA and a single pulse. 
% Initial gating variables are set to zero.

% Create vector for applied current with baseline 0.7 nA
I_app = zeros(1, length(t)) + 0.7e-9; % Baseline current

% Define stimulus parameters
pulse_start_time = 100e-3; % Pulse starts at 100 ms (converted to seconds)
pulse_duration = 5e-3; % Pulse lasts 5 ms (converted to seconds)
pulse_amplitude = 1e-9; % Pulse amplitude (1 nA)

% Calculate indices for the pulse.
pulse_start_index = floor(pulse_start_time / dt); % Start index of the pulse
pulse_end_index = pulse_start_index + floor(pulse_duration / dt); % End index of the pulse

% Apply the pulse.
I_app(pulse_start_index:pulse_end_index) = pulse_amplitude;

% Run simulation with initial gating variables set to zero.
v_init = -65e-3; % Initial membrane potential in volts (-65 mV).
m_init = 0.0;    % Initial value for 'm' (Sodium activation).
h_init = 0.0;    % Initial value for 'h' (Sodium inactivation).
n_init = 0.0;    % Initial value for 'n' (Potassium activation).

% Run the HH simulation
[v_sim, m_sim, h_sim, n_sim] = hhsim(t, I_app, v_init, m_init, h_init, n_init);

% Plot Results
figure;

% Plot Applied Current
subplot(2,1,1);
plot(t, I_app * 1e9); % Convert current to nA for plotting.
xlabel("Time (s)", 'FontSize', 12);
ylabel("I_{app} (nA)", 'FontSize', 12);
title("Part F: Applied Current with Different Initial Conditions", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Plot Membrane Potential
subplot(2,1,2);
plot(t, v_sim * 1e3); % Convert membrane potential to mV for plotting.
xlabel("Time (s)", 'FontSize', 12);
ylabel("V_{m} (mV)", 'FontSize', 12);
title("Part F: Membrane Potential with Different Initial Conditions", 'FontSize', 14, 'FontWeight', 'bold');
grid on;

% Generate plots for gating variables.
figure;
plot(t, m_sim, 'LineWidth', 1.5);
hold on;
plot(t, h_sim, '-', 'LineWidth', 1.5);
plot(t, n_sim, '--', 'LineWidth', 1.5);
xlabel("Time (s)", 'FontSize', 12);
ylabel("Gating Variables", 'FontSize', 12);
title("Gating Variables vs. Time (Part F)", 'FontSize', 14, 'FontWeight', 'bold');
legend("Sodium Activation (m)", "Sodium Inactivation (h)", "Potassium Activation (n)");
grid on;

%% Hodgkin-Huxley Simulation Function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Definitions:
% Inputs:
% t: Time vector (s)
% i_applied: Applied current (A)
% v_init: Initial membrane potential (V)
% m_init, h_init, n_init: Initial gating variables
% Outputs:
% v_sim: Simulated membrane potential (V)
% m_sim, h_sim, n_sim: Gating variables over time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [v_sim, m_sim, h_sim, n_sim] = hhsim(t, i_applied, v_init, m_init, h_init, n_init)

% Simulates the Hodgkin-Huxley model given the input time vector and applied current
global dt G_leak G_na G_k E_na E_k V_leak C_m

% Default parameters if not inputted.
if (~exist('v_init'))  % If the initial membrane potential is not provided
    v_init = -61e-3;   % Default value: -61 mV
end
if (~exist('m_init')) % If sodium activation variable is not provided
    m_init = 0;       % Default value: 0
end
if (~exist('h_init')) % If sodium inactivation variable is not provided
    h_init = 0;       % Default value: 0
end
if (~exist('n_init')) % If potassium activation variable is not provided
    n_init = 0;       % Default value: 0
end

% Setup vectors to store simulation results
v_sim = zeros(1, length(t)); % Membrane potential over time
m_sim = zeros(1, length(t)); % Sodium activation variable over time
h_sim = zeros(1, length(t)); % Sodium inactivation variable over time
n_sim = zeros(1, length(t)); % Potassium activation variable over time

% Initialize the variables at t = 0
v_sim(1) = v_init; % Initialize membrane potential
m_sim(1) = m_init; % Initialize sodium activation
h_sim(1) = h_init; % Initialize sodium inactivation
n_sim(1) = n_init; % Initialize potassium activation

% Simulation loop
for n = 1:(length(t)-1)
    % Update membrane potential
    % Calculate the ionic currents
    % Leak current: Passive flow through leak channels
    term1 = G_leak*(V_leak-v_sim(n));

    % Sodium current: Flow through voltage-gated sodium channels
    term2 = G_na*(m_sim(n)^3)*(h_sim(n))*(E_na - v_sim(n));

    % Potassium current: Flow through voltage-gated potassium channels
    term3 = G_k*(n_sim(n)^4)*(E_k - v_sim(n));

    % Use Euler's method to update the membrane potential
    v_sim(n+1) = v_sim(n) + (dt*((term1 + term2 + term3 + i_applied(n))/C_m));
    
    % Update sodium activation variable (m_sim)
    % Alpha and beta are the rate constants for m_sim
    alpha = (((10^5)*(-v_sim(n)-0.045))/(exp(100*(-v_sim(n)-0.045))-1));
    beta = (4*(10^3))*exp((-v_sim(n) - 0.07)/(0.018));
    % Rate equation for m_sim
    term1 = alpha*(1-m_sim(n));
    term2 = -beta*m_sim(n);
    m_sim(n+1) = m_sim(n) + (dt*(term1 + term2));
    
    % Update sodium inactivation variable (h_sim)
    % Alpha and beta are the rate constants for h_sim
    alpha = 70*exp(50*(-v_sim(n)-0.07));
    beta = ((10^3) / (1 + exp(100*(-v_sim(n)-0.04))));
    % Rate equation for h_sim
    term1 = alpha*(1-h_sim(n));
    term2 = -beta*h_sim(n);
    h_sim(n+1) = h_sim(n) + (dt*(term1 + term2));
    
    % Update potassium activation variable (n_sim)
    % Alpha and beta are the rate constants for n_sim
    alpha = (((10^4)*(-v_sim(n)-0.06))/(exp(100*(-v_sim(n)-0.06))-1));
    beta = 125*exp(((-v_sim(n)-0.07)/(0.08)));
    % Rate equation for n_sim
    term1 = alpha*(1-n_sim(n));
    term2 = -beta*n_sim(n);
    n_sim(n+1) = n_sim(n) + (dt*(term1 + term2)); 
end

end