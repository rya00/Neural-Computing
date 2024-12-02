% Hodgkin-Huxley Model Simulation - Part 1
clear; clc;

% Parameters
Cm = 1; % membrane capacitance (uF/cm^2)
gNa = 120; gK = 36; gL = 0.3; % Conductances (mS/cm^2)
ENa = 50; EK = -77; EL = -54.4; % Reversal potentials (mV)
V_rest = -70.2; % Resting membrane potential (mV)

% Initialize variables
V = V_rest; m = 0; h = 0; n = 0;
dt = 0.01; % Time step (ms)
tmax = 350; % Duration (ms)
time = 0:dt:tmax; % Time vector
Iapp = zeros(size(time)); % Applied current (nA)

% Allocate arrays
Vm = zeros(size(time)); m_trace = Vm; h_trace = Vm; n_trace = Vm;

% Gating variable functions
alpha_m = @(V) 0.1*(V+40)/(1 - exp(-(V+40)/10));
beta_m = @(V) 4*exp(-0.0556*(V+65));
alpha_h = @(V) 0.07*exp(-0.05*(V+65));
beta_h = @(V) 1/(1 + exp(-0.1*(V+35)));
alpha_n = @(V) 0.01*(V+55)/(1 - exp(-(V+55)/10));
beta_n = @(V) 0.125*exp(-0.0125*(V+65));

% Simulation loop
for i = 1:length(time)
    % Ionic currents
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    IL = gL * (V - EL);

    % Update membrane potential
    dVdt = (Iapp(i) - INa - IK - IL) / Cm;
    V = V + dt * dVdt;

    % Update gating variables
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    m = m + dt * dm;
    h = h + dt * dh;
    n = n + dt * dn;

    % Store values
    Vm(i) = V; m_trace(i) = m; h_trace(i) = h; n_trace(i) = n;
end

% Plot results
figure;
plot(time, Vm);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Hodgkin-Huxley Model: Resting Potential');

%% Question 2
% Add step current
Iapp(time > 100 & time <= 200) = 0.22; % 0.22 nA step current

% Simulation loop
for i = 1:length(time)
    % Ionic currents
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    IL = gL * (V - EL);

    % Update membrane potential
    dVdt = (Iapp(i) - INa - IK - IL) / Cm;
    V = V + dt * dVdt;

    % Update gating variables
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    m = m + dt * dm;
    h = h + dt * dh;
    n = n + dt * dn;

    % Store values
    Vm(i) = V; m_trace(i) = m; h_trace(i) = h; n_trace(i) = n;
end

% Plot applied current and membrane potential
figure;
subplot(2, 1, 1);
plot(time, Iapp);
xlabel('Time (ms)');
ylabel('Applied Current (nA)');
title('Step Current Applied');

subplot(2, 1, 2);
plot(time, Vm);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential Response to Step Current');

%% Question 3

% Generate 10 pulses of 5 ms each with varying delays
pulse_duration = 5; % ms
delay = 15; % Set delay between pulses (adjust as needed)
Iapp(:) = 0; % Reset applied current

for j = 0:9
    start_idx = find(time >= (100 + j * delay), 1);
    end_idx = find(time >= (100 + j * delay + pulse_duration), 1);
    Iapp(start_idx:end_idx) = 0.22; % Pulse amplitude
end

% Simulation loop
for i = 1:length(time)
    % Ionic currents
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    IL = gL * (V - EL);

    % Update membrane potential
    dVdt = (Iapp(i) - INa - IK - IL) / Cm;
    V = V + dt * dVdt;

    % Update gating variables
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    m = m + dt * dm;
    h = h + dt * dh;
    n = n + dt * dn;

    % Store values
    Vm(i) = V; m_trace(i) = m; h_trace(i) = h; n_trace(i) = n;
end

% Plot applied current and membrane potential
figure;
subplot(2, 1, 1);
plot(time, Iapp);
xlabel('Time (ms)');
ylabel('Applied Current (nA)');
title('Pulsed Current Applied');

subplot(2, 1, 2);
plot(time, Vm);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Membrane Potential Response to Pulsed Current');

%% Question 4

% Baseline current of 0.6 nA with inhibitory pulses
Iapp(:) = 0.6; % Baseline current
for j = 0:9
    start_idx = find(time >= (100 + j * 20), 1);
    end_idx = find(time >= (100 + j * 20 + 5), 1);
    Iapp(start_idx:end_idx) = 0; % Bring current to zero
end

% Set initial conditions
V = -65; m = 0.05; h = 0.5; n = 0.35;

% Simulation loop
for i = 1:length(time)
    % Ionic currents
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    IL = gL * (V - EL);

    % Update membrane potential
    dVdt = (Iapp(i) - INa - IK - IL) / Cm;
    V = V + dt * dVdt;

    % Update gating variables
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    m = m + dt * dm;
    h = h + dt * dh;
    n = n + dt * dn;

    % Store values
    Vm(i) = V; m_trace(i) = m; h_trace(i) = h; n_trace(i) = n;
end

% Plot applied current and membrane potential
figure;
subplot(2, 1, 1);
plot(time, Iapp);
xlabel('Time (ms)');
ylabel('Applied Current (nA)');
title('Inhibitory Pulses Applied');

subplot(2, 1, 2);
plot(time, Vm);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Response to Inhibitory Pulses');

%% Question 5

% Part 5: Baseline 0.65 nA with 1 nA excitatory pulse
Iapp(:) = 0.65; % Baseline current
start_idx = find(time >= 100, 1);
end_idx = find(time >= 105, 1);
Iapp(start_idx:end_idx) = 1.65; % 1 nA pulse added to baseline

% Part 6: Baseline 0.7 nA with 1 nA pulse and different initial conditions
% Reset initial conditions
V = -65; 
m = 0; 
h = 0; 
n = 0;
Iapp(:) = 0.7; % Baseline current
Iapp(start_idx:end_idx) = 1.7; % 1 nA pulse added to baseline

for i = 1:length(time)
    % Ionic currents
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    IL = gL * (V - EL);

    % Update membrane potential
    dVdt = (Iapp(i) - INa - IK - IL) / Cm;
    V = V + dt * dVdt;

    % Update gating variables
    dm = alpha_m(V) * (1 - m) - beta_m(V) * m;
    dh = alpha_h(V) * (1 - h) - beta_h(V) * h;
    dn = alpha_n(V) * (1 - n) - beta_n(V) * n;
    m = m + dt * dm;
    h = h + dt * dh;
    n = n + dt * dn;

    % Store values
    Vm(i) = V; m_trace(i) = m; h_trace(i) = h; n_trace(i) = n;
end

% Plot applied current and membrane potential
figure;
subplot(2, 1, 1);
plot(time, Iapp);
xlabel('Time (ms)');
ylabel('Applied Current (nA)');
title('Excitatory Pulse Applied');

subplot(2, 1, 2);
plot(time, Vm);
xlabel('Time (ms)');
ylabel('Membrane Potential (mV)');
title('Response to Excitatory Pulse');

