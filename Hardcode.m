% Hodgkin-Huxley Neuron Model Simulation
clear all; close all;

% Parameters
dt = 2e-8;          % Time step (s)
tmax = 0.35;        % Simulation duration (s)
t = 0:dt:tmax;      % Time vector

% Constants (in SI units)
V_L = -0.060;       % Leak reversal potential (V)
E_Na = 0.045;       % Sodium reversal potential (V)
E_K = -0.082;       % Potassium reversal potential (V)
G_L = 30e-9;        % Leak conductance (S)
G_Na = 12e-6;       % Sodium conductance (S)
G_K = 3.6e-6;       % Potassium conductance (S)
Cm = 100e-12;       % Membrane capacitance (F)

% Default initial conditions
V0 = -0.065;        % Initial membrane potential (V)
m0 = 0.05;          % Initial sodium activation gating variable
h0 = 0.5;           % Initial sodium inactivation gating variable
n0 = 0.35;          % Initial potassium activation gating variable

% Applied current parameters
istart = 100e-3;    % Start time of applied current (s)
ilength = 5e-3;     % Duration of applied current pulse (s)
Ibase = 0;          % Baseline current (A)

% Initialize the figures
figure;

% --- Part 2: Step current ---
% Set parameters specific to Part 2
ilength = 100e-3;   % 100 ms step current
Ie = 0.22e-9;       % 0.22 nA step amplitude

% Applied current for Part 2
Iapp_part2 = Ibase * ones(size(t));
Iapp_part2(t > 100e-3 & t <= 200e-3) = Ie;

% Run simulation for Part 2
[V_part2, Iapp_part2_plot] = simulate_HH(Iapp_part2, V0, m0, h0, n0);

% --- Part 3: Pulses ---
% Set parameters specific to Part 3
Npulses = 10;       % 10 pulses
Ie = 0.22e-9;       % 0.22 nA pulse amplitude
pulsesep = 18e-3;   % Pulse separation

% Applied current for Part 3
Iapp_part3 = Ibase * ones(size(t));
for pulse = 1:Npulses
    start_idx = round((istart + (pulse - 1) * pulsesep) / dt) + 1;
    stop_idx = round((istart + (pulse - 1) * pulsesep + ilength) / dt) + 1;
    Iapp_part3(start_idx:stop_idx) = Ie;
end

% Run simulation for Part 3
[V_part3, Iapp_part3_plot] = simulate_HH(Iapp_part3, V0, m0, h0, n0);

% --- Part 4: Inhibitory pulses ---
% Set parameters specific to Part 4
Npulses = 10;       % 10 inhibitory pulses
Ibase = 0.6e-9;     % Baseline current
Ie = 0;             % Bring current to zero during pulse

% Applied current for Part 4
Iapp_part4 = Ibase * ones(size(t));
for pulse = 1:Npulses
    start_idx = round((istart + (pulse - 1) * pulsesep) / dt) + 1;
    stop_idx = round((istart + (pulse - 1) * pulsesep + ilength) / dt) + 1;
    Iapp_part4(start_idx:stop_idx) = Ie;
end

% Run simulation for Part 4
[V_part4, Iapp_part4_plot] = simulate_HH(Iapp_part4, V0, m0, h0, n0);

% --- Part 5: Excitatory current ---
% Set parameters specific to Part 5
Ibase = 0.65e-9;    % Baseline current
Ie = 1e-9;          % 1 nA excitatory pulse

% Applied current for Part 5
Iapp_part5 = Ibase * ones(size(t));
start_idx = round(100e-3 / dt) + 1; % 100 ms start time
stop_idx = round((100e-3 + ilength) / dt) + 1; % 5 ms duration pulse
Iapp_part5(start_idx:stop_idx) = Ie;

% Run simulation for Part 5
[V_part5, Iapp_part5_plot] = simulate_HH(Iapp_part5, V0, m0, h0, n0);

% --- Part 6: Excitatory current with zero initial conditions ---
% Set parameters specific to Part 6
Ibase = 0.65e-9;    % Baseline current
Ie = 1e-9;          % 1 nA excitatory pulse
m0 = 0;             % Zero initial sodium activation
h0 = 0;             % Zero initial sodium inactivation
n0 = 0;             % Zero initial potassium activation

% Applied current for Part 6
Iapp_part6 = Ibase * ones(size(t));
start_idx = round(100e-3 / dt) + 1; % 100 ms start time
stop_idx = round((100e-3 + ilength) / dt) + 1; % 5 ms duration pulse
Iapp_part6(start_idx:stop_idx) = Ie;

% Run simulation for Part 6
[V_part6, Iapp_part6_plot] = simulate_HH(Iapp_part6, V0, m0, h0, n0);

%% Plot results
subplot(5, 2, 1);
plot(t, Iapp_part2 * 1e9); % Applied current (nA)
ylabel('I_{app} (nA)');
title('Part 2: Step Current');
subplot(5, 2, 2);
plot(t, V_part2 * 1e3); % Membrane potential (mV)
ylabel('V_m (mV)');
xlabel('Time (s)');

subplot(5, 2, 3);
plot(t, Iapp_part3 * 1e9); % Applied current (nA)
ylabel('I_{app} (nA)');
title('Part 3: Pulsed Current');
subplot(5, 2, 4);
plot(t, V_part3 * 1e3); % Membrane potential (mV)
ylabel('V_m (mV)');
xlabel('Time (s)');

subplot(5, 2, 5);
plot(t, Iapp_part4 * 1e9); % Applied current (nA)
ylabel('I_{app} (nA)');
title('Part 4: Inhibitory Pulses');
subplot(5, 2, 6);
plot(t, V_part4 * 1e3); % Membrane potential (mV)
ylabel('V_m (mV)');
xlabel('Time (s)');

subplot(5, 2, 7);
plot(t, Iapp_part5 * 1e9); % Applied current (nA)
ylabel('I_{app} (nA)');
title('Part 5: Excitatory Current');
subplot(5, 2, 8);
plot(t, V_part5 * 1e3); % Membrane potential (mV)
ylabel('V_m (mV)');
xlabel('Time (s)');

subplot(5, 2, 9);
plot(t, Iapp_part6 * 1e9); % Applied current (nA)
ylabel('I_{app} (nA)');
title('Part 6: Excitatory Current (Zero Initial Conditions)');
subplot(5, 2, 10);
plot(t, V_part6 * 1e3); % Membrane potential (mV)
ylabel('V_m (mV)');
xlabel('Time (s)');

%% Function to simulate the Hodgkin-Huxley model
function [V, Iapp] = simulate_HH(Iapp, V0, m0, h0, n0)
    dt = 2e-8;          % Time step (s)
    tmax = 0.35;        % Simulation duration (s)
    t = 0:dt:tmax;      % Time vector

    % Constants (in SI units)
    V_L = -0.060;       % Leak reversal potential (V)
    E_Na = 0.045;       % Sodium reversal potential (V)
    E_K = -0.082;       % Potassium reversal potential (V)
    G_L = 30e-9;        % Leak conductance (S)
    G_Na = 12e-6;       % Sodium conductance (S)
    G_K = 3.6e-6;       % Potassium conductance (S)
    Cm = 100e-12;       % Membrane capacitance (F)

    % Initialize variables
    V = zeros(size(t)); V(1) = V0;
    m = zeros(size(t)); m(1) = m0;
    h = zeros(size(t)); h(1) = h0;
    n = zeros(size(t)); n(1) = n0;
    I_Na = zeros(size(t)); I_K = zeros(size(t)); I_L = zeros(size(t)); Itot = zeros(size(t));

    % Simulation loop
    for i = 2:length(t)
        Vm = V(i-1);  % Membrane potential

        % Gating variable rates
        alpha_m = (1e5 * (-Vm - 0.045)) / (exp(100 * (-Vm - 0.045)) - 1);
        beta_m = 4000 * exp((-Vm - 0.070) / 0.018);
        alpha_h = 70 * exp(50 * (-Vm - 0.070));
        beta_h = 1000 / (1 + exp(100 * (-Vm - 0.040)));
        alpha_n = (1e4 * (-Vm - 0.060)) / (exp(100 * (-Vm - 0.060)) - 1);
        beta_n = 125 * exp((-Vm - 0.070) / 0.08);

        % Time constants and steady-state values
        tau_m = 1 / (alpha_m + beta_m); m_inf = alpha_m / (alpha_m + beta_m);
        tau_h = 1 / (alpha_h + beta_h); h_inf = alpha_h / (alpha_h + beta_h);
        tau_n = 1 / (alpha_n + beta_n); n_inf = alpha_n / (alpha_n + beta_n);

        % Update gating variables
        m(i) = m(i-1) + (m_inf - m(i-1)) * dt / tau_m;
        h(i) = h(i-1) + (h_inf - h(i-1)) * dt / tau_h;
        n(i) = n(i-1) + (n_inf - n(i-1)) * dt / tau_n;

        % Ionic currents
        I_Na(i) = G_Na * m(i)^3 * h(i) * (E_Na - Vm);
        I_K(i) = G_K * n(i)^4 * (E_K - Vm);
        I_L(i) = G_L * (V_L - Vm);

        % Total current and membrane potential update
        Itot(i) = I_Na(i) + I_K(i) + I_L(i) + Iapp(i);
        V(i) = Vm + (Itot(i) * dt / Cm);
    end
end
