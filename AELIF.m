% Adaptive Exponential LIF (AELIF) Model Parameters
E_L = -75e-3;          % Leak reversal potential (V)
V_th = -50e-3;         % Threshold potential (V)
V_reset = -80e-3;      % Reset potential (V)
Delta_th = 2e-3;       % Exponential threshold factor (V)
g_L = 10e-9;           % Leak conductance (S)
C_m = 100e-12;         % Membrane capacitance (F)
a = 2e-9;              % Subthreshold adaptation (S)
b = 20e-12;            % Spike-triggered adaptation increment (A)
tau_SRA = 0.2;         % Adaptation time constant (s)

% Simulation Parameters and Time Vector
dt = 0.0001;           % Time step (s)
T = 1.5;               % Total simulation time (s)
time = 0:dt:T;         % Time vector

% Define input current pulse (500 pA from 0.5s to 1.0s)
I_app = 500e-12;       % Applied current (A)
I_input = zeros(size(time));
I_input(time >= 0.5 & time <= 1.0) = I_app;

% Initialize Variables
V = E_L;               % Initial Membrane potential (V)
I_SRA = 0;             % Initial Adaptation current (A)
V_trace = zeros(size(time)); % Store Membrane potential
I_SRA_trace = zeros(size(time)); % Store Adaptation conductance

% Simulation Loop
for t = 1:length(time)
    % Calculate Currents
    I_leak = g_L * (E_L - V); % Leak current (A)
    I_exp = g_L * Delta_th * exp((V - V_th) / Delta_th); % Exponential current (A)

    % Update Membrane Potential Using Euler Method
    dVdt = (I_leak + I_exp + I_input(t) - I_SRA) / C_m;
    V = V + dt * dVdt;

    % Check for Spike and Reset Threshold if Reached
    if V > V_th
        V = V_reset;
        I_SRA = I_SRA + b; % Spike-triggered adaptation
    end

    % Update Adaptation Current
    dI_SRA_dt = (a * (V - E_L) - I_SRA) / tau_SRA;
    I_SRA = I_SRA + dt * dI_SRA_dt;

    % Store Results
    V_trace(t) = V;
    I_SRA_trace(t) = I_SRA;
end

% Plot Results for Part A
figure;
subplot(2, 1, 1);
plot(time, I_input * 1e12); % Convert to pA for plotting
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Input Current (pA)', 'FontWeight', 'bold');
title('Input Current vs Time', 'FontWeight', 'bold');
grid on;

subplot(2, 1, 2);
plot(time, V_trace * 1e3); % Convert to mV for plotting
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Membrane Potential (mV)', 'FontWeight', 'bold');
title('Membrane Potential vs Time', 'FontWeight', 'bold');
grid on;

sgtitle('AELIF Model Simulation: Part A', 'FontSize', 16, 'FontWeight', 'bold');

% Parameters for Part B
T_B = 5;               % Total simulation time (s)
time_B = 0:dt:T_B;     % Time vector
I_range = linspace(0, 1e-9, 20); % Range of applied currents (A)

% Storage for Results
initial_firing_rates = zeros(size(I_range)); % Initial firing rates (Hz)
steady_firing_rates = zeros(size(I_range));  % Steady-state firing rates (Hz)

% Simulation for f-I Curve
V_B = E_L * ones(length(time_B), length(I_range)); % Membrane potentials
I_SRA_B = zeros(length(time_B), length(I_range));  % Adaptation currents
spike_times_B = cell(1, length(I_range));

for t_B = 2:length(time_B)
    % Update Membrane Potential
    I_leak_B = g_L * (E_L - V_B(t_B-1, :)); % Leak current (A)
    I_exp_B = g_L * Delta_th * exp((V_B(t_B-1, :) - V_th) / Delta_th); % Exponential current (A)
    dVdt_B = (I_leak_B + I_exp_B + I_range - I_SRA_B(t_B-1, :)) / C_m;
    V_B(t_B, :) = V_B(t_B-1, :) + dt * dVdt_B;

    % Check for Spikes
    spike_indices = V_B(t_B, :) > V_th;
    V_B(t_B, spike_indices) = V_reset; % Reset potential
    I_SRA_B(t_B, spike_indices) = I_SRA_B(t_B-1, spike_indices) + b;

    % Record Spike Times
    for i = find(spike_indices)
        spike_times_B{i}(end + 1) = t_B * dt;
    end

    % Update Adaptation Current
    dI_SRA_dt_B = (a * (V_B(t_B, :) - E_L) - I_SRA_B(t_B-1, :)) / tau_SRA;
    I_SRA_B(t_B, :) = I_SRA_B(t_B-1, :) + dt * dI_SRA_dt_B;
end

% Calculate Firing Rates
for i = 1:length(I_range)
    if length(spike_times_B{i}) > 1
        ISIs_B = diff(spike_times_B{i}); % Interspike intervals (s)
        initial_firing_rates(i) = 1 / ISIs_B(1); % Initial firing rate (Hz)
        steady_firing_rates(i) = 1 / mean(ISIs_B(end - min(10, length(ISIs_B)-1):end)); % Steady-state rate (Hz)
    else
        initial_firing_rates(i) = 0;
        steady_firing_rates(i) = 0;
    end
end

% Plot f-I Curve for Part B
figure;
hold on;
plot(I_range * 1e12, steady_firing_rates, 'k-', 'DisplayName', 'Steady-State Firing Rate'); % pA
plot(I_range * 1e12, initial_firing_rates, 'ko', 'DisplayName', 'Initial Firing Rate');
xlabel('Applied Current (pA)', 'FontWeight', 'normal');
ylabel('Firing Rate (Hz)', 'FontWeight', 'normal');
title('f-I Curve for AELIF Model', 'FontWeight', 'normal');
legend('show');
grid on;
hold off;
