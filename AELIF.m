% Adaptive Exponential LIF (AELIF) Model Parameters
E_L = -75; % Leak reversal potential (mV)
V_th = -50; % Threshold potential (mV)
V_reset = -80; % Reset potential (mV)
Delta_th = 2; % Exponential threshold factor (mV)
g_L = 10; % Leak conductance (nS)
C_m = 100; % Membrane capacitance (pF)
a = 2; % Subthreshold adaptation (nS)
b = 0.02; % Spike-triggered adaptation increment (nA)
tau_SRA = 200; % Adaptation time constant (ms)

% Simulation Parameters and Create Time Vector
dt = 0.1; % Time step (ms)
T = 1500; % Total simulation time (ms)
time = 0:dt:T; % Time vector

% Define input current pulse (500 pA from 0.5s to 1.0s)
I_app = 500; % Applied current (pA)
I_input = zeros(size(time));
I_input(time >= 500 & time <= 1000) = I_app;

% Initial variables
V = E_L; % Initial Membrane potential (mV)
I_SRA = 0; % Initial Adaptation current (nA)
V_trace = zeros(size(time)); % Store Membrane potential
I_SRA_trace = zeros(size(time)); % Store Adaptation conductance

% Simulation loop
for t = 1:length(time)
    % Calculate currents
    I_leak = g_L * (E_L - V); % Leak current (pA)
    I_exp = g_L * Delta_th * exp((V - V_th) / Delta_th); % Exponential current (pA)

    % Update membrane potential using Euler method
    dVdt = (I_leak + I_exp + I_input(t) - I_SRA) / C_m;
    V = V + dt * dVdt;

    % Check for spike and reset threshold if reached
    if V > V_th
        V = V_reset;
        I_SRA = I_SRA + b; % Spike-triggered adaptation
    end

    % Update adaptation current
    dI_SRA_dt = (a * (V - E_L) - I_SRA) / tau_SRA;
    I_SRA = I_SRA + dt * dI_SRA_dt;

    % Store results
    V_trace(t) = V;
    I_SRA_trace(t) = I_SRA;
end

% Set consistent font size and style for all plots
set(0, 'DefaultAxesFontSize', 12)
set(0, 'DefaultAxesFontName', 'Arial')

% Plot results for Part A
figure;
subplot(2, 1, 1);
plot(time, I_input);
xlabel('Time (ms)', 'FontWeight', 'bold');
ylabel('Input Current (pA)', 'FontWeight', 'bold');
title('Input Current vs Time', 'FontWeight', 'bold');
grid on;

subplot(2, 1, 2);
plot(time, V_trace);
xlabel('Time (ms)', 'FontWeight', 'bold');
ylabel('Membrane Potential (mV)', 'FontWeight', 'bold');
title('Membrane Potential vs Time', 'FontWeight', 'bold');
grid on;

sgtitle('AELIF Model Simulation', 'FontSize', 16, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 800, 600]);

% Parameters for Part B
dt_B = 0.1; % Time step (ms)
T_B = 5000; % Total simulation time (ms)
time_B = 0:dt_B:T_B; % Time vector

I_range = linspace(0, 1000, 20); % Range of applied currents (pA)
firing_rates = zeros(size(I_range)); % Steady-state firing rates
initial_firing_rates = zeros(size(I_range)); % Initial firing rates

% Vectorized simulation for f-I curve
V_B = E_L * ones(length(time_B), length(I_range));
I_SRA_B = zeros(length(time_B), length(I_range));
spike_times_B = cell(1, length(I_range));

for t_B = 2:length(time_B)
    % Update membrane potential
    I_leak_B = g_L * (E_L - V_B(t_B-1, :));
    I_exp_B = g_L * Delta_th * exp((V_B(t_B-1, :) - V_th) / Delta_th);
    dVdt_B = (I_leak_B + I_exp_B + I_range - I_SRA_B(t_B-1, :)) / C_m;
    V_B(t_B, :) = V_B(t_B-1, :) + dt_B * dVdt_B;
    
    % Check for spikes
    spike_indices = V_B(t_B, :) > V_th;
    V_B(t_B, spike_indices) = V_reset;
    I_SRA_B(t_B, spike_indices) = I_SRA_B(t_B-1, spike_indices) + b;
    
    % Record spike times
    for i = find(spike_indices)
        spike_times_B{i}(end + 1) = t_B * dt_B;
    end
    
    % Update adaptation current
    dI_SRA_dt_B = (a * (V_B(t_B, :) - E_L) - I_SRA_B(t_B-1, :)) / tau_SRA;
    I_SRA_B(t_B, :) = I_SRA_B(t_B-1, :) + dt_B * dI_SRA_dt_B;
end

% Calculate firing rates
for i = 1:length(I_range)
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

% Plot the f-I curve
figure;
plot(I_range, firing_rates, '-o', 'DisplayName', 'Steady-State Firing Rate');
hold on;
plot(I_range, initial_firing_rates, 'x', 'DisplayName', 'Initial Firing Rate');
xlabel('Applied Current (pA)');
ylabel('Firing Rate (Hz)');
title('f-I Curve for AELIF Model');
legend;
grid on;

% Check if firing rate range is adequate
max_firing_rate = max(firing_rates);
if max_firing_rate < 50
    warning('Maximum firing rate (%.2f Hz) is below 50 Hz. Consider increasing the current range.', max_firing_rate);
elseif max_firing_rate > 60
    warning('Maximum firing rate (%.2f Hz) exceeds 60 Hz. Consider decreasing the current range.', max_firing_rate);
else
    disp('Firing rate range is appropriate (0 to 50+ Hz).');
end

% Comment on results
disp('The initial firing rates are generally higher than the steady-state firing rates due to adaptation.');
disp('At higher applied currents, the difference between initial and steady-state rates increases, reflecting stronger adaptation effects.');
