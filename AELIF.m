%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Exponential Leaky Integrate-and-Fire (AELIF) Neuron Model
% This script simulates the AELIF neuron and analyzes its behavior:
% Part A: Simulates the effect of a 500 pA pulse on the neuron for 1.5s.
% Part B: Examines firing behavior across a range of applied currents and computes the f-I curve.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization for Part A 
clear all; % Clear all variables from the workspace
close all; % Close all open figures

% Adaptive Exponential LIF (AELIF) Model Parameters
E_l = -75e-3;          % Leak reversal potential (V)
V_th = -50e-3;         % Threshold potential (V)
V_reset = -80e-3;      % Reset potential (V)
Delta_th = 2e-3;       % Exponential threshold factor (V)
G_l = 10e-9;           % Leak conductance (S)
C_m = 100e-12;         % Membrane capacitance (F)
a = 2e-9;              % Subthreshold adaptation (S)
b = 20e-12;            % Spike-triggered adaptation increment (A)
Tau_SRA = 0.2;         % Adaptation time constant (s)

%% Simulation Parameters and Time Vector for Part A
dt = 0.0001;           % Time step (s)
T = 1.5;               % Total simulation time (s)
time = 0:dt:T;         % Time vector

% Define input current pulse (500 pA from 0.5s to 1.0s)
I_input = zeros(size(time)); % Initialize applied current as zero
I_input(time >= 0.5 & time <= 1.0) = 500e-12; % Apply 500 pA during the specified interval

% Initialize Variables
V = E_l;               % Initial Membrane potential (V)
I_SRA = 0;             % Initial Adaptation current (A)
V_trace = zeros(size(time)); % Store Membrane potential
I_SRA_trace = zeros(size(time)); % Store Adaptation conductance

% Simulation Loop
for t = 1:length(time)
    % Calculate Currents
    I_leak = G_l * (E_l - V); % Leak current: passive flow towards leak potential
    I_exp = G_l * Delta_th * exp((V - V_th) / Delta_th); % Exponential current

    % Update Membrane Potential Using Euler Method
    dVdt = (I_leak + I_exp + I_input(t) - I_SRA) / C_m; % Total effect on membrane potential
    V = V + dt * dVdt; % Euler integration step

    % Check for Spike and Reset Threshold if Reached
    if V > V_th
        V = V_reset; % Reset membrane potential after spike
        I_SRA = I_SRA + b; % Spike-triggered adaptation
    end

    % Update Adaptation Current
    dI_SRA_dt = (a * (V - E_l) - I_SRA) / Tau_SRA; % Adaptation dynamics
    I_SRA = I_SRA + dt * dI_SRA_dt; % Euler integration step for adaptation current

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

% Simulation Parameters
%% Part B: f-I Curve Simulation
I_base = 0e-9;              % Baseline current
V_max = 50e-3;              % Level of voltage to detect a spike

% Simulation Parameters and Time Vector for Part B
dt = 2e-6;                   % Time step (s)
t_max = 5;                   % Total simulation time (s)
t_vector = 0:dt:t_max;       % Time vector

I = I_base*ones(size(t_vector));          % Initialize input current vector
current_range = [0.15:0.005:0.5]*1e-9;    % Range of applied currents

% Initialize Storage for Results
initial_rate = zeros(size(current_range)); % Store initial spike rate
final_rate = zeros(size(current_range));   % Store steady-state spike rate
single_spike = zeros(size(current_range)); % Track trials with a single spike
mean_V = zeros(size(current_range));       % Average membrane potential

% Simulation Loop for Part B
for i = 1:length(current_range)
    I_app_b = current_range(i); % Current level being tested
    I(:) = I_app_b; % Apply constant current across simulation
    
    % Initialize variables
    v = zeros(size(t_vector)); % Membrane potential
    v(1) = E_l; % Set value of initial membrane potential
    I_sra = zeros(size(t_vector)); % Initialize adaptation current
    spikes = zeros(size(t_vector)); % Initialize vector to store spikes
    
    for j = 1:length(t_vector)-1        
        if ( v(j) > V_max )
            v(j) = V_reset; % Reset membrane potential
            I_sra(j) = I_sra(j) + b; % Increase the adaptation current by b
            spikes(j) = 1; % Record spike
        end
        
        % Update Membrane Potential
        v(j+1) = v(j) + dt*( G_l*(E_l-v(j) + Delta_th*exp((v(j)-V_th)/Delta_th) ) - I_sra(j) + I(j))/C_m;
        
        % Update Adaptation Current
        I_sra(j+1) = I_sra(j) + dt*( a*(v(j)-E_l) - I_sra(j) )/Tau_SRA;
        
    end
    
    % Analyze Spike Data
    spike_times = dt*find(spikes);        % Extract the spike times in s
    if ( length(spike_times) > 1 )           
        ISIs = diff(spike_times);         % Interspike intervals
        initial_rate(i) = 1/ISIs(1);      % Initial spike rate
        if ( length(ISIs) > 1 )             
            final_rate(i) = 1/ISIs(end);  % Steady-state spike rate
        end
    else
        if ( length(spike_times) == 1 )      
            single_spike(i) = 1;          % Record single spike occurrence
        end
    end
    
    mean_V(i) = mean(v); % Average membrane potential
end

% Plot f-I Curve for Part B
figure;
hold on;                            % Allow many plots on same graph
plot(1e9*current_range, final_rate, 'b', 'LineWidth', 2);     
ISIindices = find(initial_rate);     % Find points where a first ISI exists
plot(1e9*current_range(ISIindices), initial_rate(ISIindices), 'ro', 'MarkerSize', 6);
ISIindices = find(single_spike);     % Find points with just one spike (no ISI)
plot(1e9*current_range(ISIindices), 0*single_spike(ISIindices), '*k', 'MarkerSize', 8);
xlabel('Iapp (nA)');                % Label x-axis
ylabel('Spike Rate (Hz)');          % Label y-axis
legend('Final Rate', '1/ISI(1)', 'Single spike');
hold off;
grid on;

sgtitle('AELIF f-I Curve (Part B)', 'FontSize', 16, 'FontWeight', 'bold');