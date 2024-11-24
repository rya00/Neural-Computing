% Code for adaptive LIF neuron model
clear all;

% Parameters of the LIF model
E_L = -75; % Leak reversal potential (mV)
V_th = -50; % Threshold potential (mV)
V_reset = -80; % Reset potential (mV)
R_m = 100; % Membrane resistance (MΩ)
C_m = 100; % Membrane capacitance (pF)
E_k = -80; % Adaptation reversal potential (mV)
Delta_GSRA = 1; % Spike-triggered adaptation increment (nS)
tau_SRA = 200; % Adaptation time constant (ms)
I_app = 500; % Applied current (pA)

% Simulation Parameters and Create Time Vector for Mode Neuron
dt = 0.1; % Time step (ms)
T = 1500; % Total simulation time (ms) (1.5s)
time = 0:dt:T; % Time vector

% Generate input current pulse (500 pA from 500 ms to 1000 ms)
I_input = zeros(size(time));
I_input(time >= 500 & time <= 1000) = I_app; 

% Initialise variables
V = E_L; % Initial Membrane potential (mV)
GSRA = 0; % Initial Adaptation conductance (nS)
V_trace = zeros(size(time)); % Store Membrane potential
GSRA_trace = zeros(size(time)); % Store Adaptation conductance

% Simulation loop
for t = 1:length(time)
    % Calculate current according to given formula
    I_leak = (E_L - V) / R_m; % Leak current (pA) calculation
    I_adapt = GSRA * (E_k - V); % Adaptation current (pA)

    % Update membrane potential
    dVdt = (I_leak + I_adapt + I_input(t)) / C_m; % Membrane potential derivative (mV/ms)
    V = V + dt * dVdt; 

    % Check for spike and reset if threshold is reached
    if V > V_th
        V = V_reset; 
        GSRA = GSRA + Delta_GSRA; 
    end

    % Update adaptation conductance
    dGSRA_dt = -GSRA / tau_SRA; % Adaptation conductance decay (nS/ms)
    GSRA = GSRA + dt * dGSRA_dt; % Update adaptation conductance

    % Store results
    V_trace(t) = V;
    GSRA_trace(t) = GSRA;
end

% Set consistent font size and style for all plots
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultAxesFontName', 'Arial')

% Plot results for Part A
figure;
subplot(3, 1, 1);
plot(time, I_input);
xlabel('Time (ms)', 'FontWeight', 'bold');
ylabel('Input Current (pA)', 'FontWeight', 'bold');
title('Input Current vs Time', 'FontWeight', 'bold');
grid on;

subplot(3, 1, 2);
plot(time, V_trace);
xlabel('Time (ms)', 'FontWeight', 'bold');
ylabel('Membrane Potential (mV)', 'FontWeight', 'bold');
title('Membrane Potential vs Time', 'FontWeight', 'bold');
grid on;

subplot(3, 1, 3);
plot(time, GSRA_trace);
xlabel('Time (ms)', 'FontWeight', 'bold');
ylabel('Adaptation Conductance (nS)', 'FontWeight', 'bold');
title('Adaptation Conductance vs Time', 'FontWeight', 'bold');
grid on;

sgtitle('LIF Model with Adaptation Current', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 800, 600]);

% Simulation parameters for model with 20 different levels of constant applied
dt_B = 0.1;     % Time step (ms)
T_B = 5000;     % Total simulation time (ms)
time_B = 0:dt_B:T_B; % Time vector

% Adjust range of applied currents
I_range = linspace(0, 1500, 20); 

% Storage for firing rates and ISIs
firing_rates = zeros(size(I_range));
initial_firing_rates = zeros(size(I_range));

% Vectorized simulation for f-I curve
V_B = E_L * ones(length(time_B), 1);
GSRA_B = zeros(length(time_B), 1);
spike_times_B = cell(1, length(I_range));

for i = 1:length(I_range)
    I_app_B = I_range(i);
    
    for t_B = 2:length(time_B)
        % Update membrane potential
        I_leak_B = (E_L - V_B(t_B-1)) / R_m;
        I_adapt_B = GSRA_B(t_B-1) * (E_k - V_B(t_B-1));
        dVdt_B = (I_leak_B + I_adapt_B + I_app_B) / C_m;
        V_B(t_B) = V_B(t_B-1) + dt_B * dVdt_B;

        % Check for spike
        if V_B(t_B) > V_th
            spike_times_B{i}(end + 1) = t_B * dt_B;
            V_B(t_B) = V_reset;
            GSRA_B(t_B) = GSRA_B(t_B-1) + Delta_GSRA;
        else
            GSRA_B(t_B) = GSRA_B(t_B-1) + dt_B * (-GSRA_B(t_B-1) / tau_SRA);
        end
    end
    
    % Calculate firing rates
    if length(spike_times_B{i}) > 1
        ISIs_B = diff(spike_times_B{i});
        initial_firing_rates(i) = 1000 / ISIs_B(1);
        steady_state_ISIs = ISIs_B(end - min(10, length(ISIs_B)-1):end);
        firing_rates(i) = 1000 / mean(steady_state_ISIs);
    else
        initial_firing_rates(i) = 0;
        firing_rates(i) = 0;
    end
end

% Plot f-I Curve
figure;
plot(I_range, firing_rates, '-o', 'DisplayName', 'Steady-State Firing Rate'); % Continuous curve
hold on;
plot(I_range, initial_firing_rates, 'x', 'DisplayName', 'Initial Firing Rate'); % Individual points
xlabel('Applied Current (pA)');
ylabel('Firing Rate (Hz)');
title('f-I Curve for LIF Model with Adaptation');
legend;
grid on;

% Check if firing rate range is adequate
max_firing_rate = max(firing_rates);
if max_firing_rate < 50
    warning('Maximum firing rate is below 50 Hz. Consider increasing the current range.', max_firing_rate);
elseif max_firing_rate > 60
    warning('Maximum firing rate exceeds 60 Hz. Consider decreasing the current range.', max_firing_rate);
else
    disp('Firing rate range is appropriate (0 to 50+ Hz).');
end

