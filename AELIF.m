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

% Simulation Parameters
%% Part B: f-I Curve Simulation
I_base = 0e-9;                 % Baseline current (A)

V_max = 50e-3;              % level of voltage to detect a spike

dt = 2e-6;                  % time-step in sec
t_max = 5;                   % maximum time in sec
t_vector = 0:dt:t_max;        % vector of all the time points

I = I_base*ones(size(t_vector)); % baseline current added to all time points

current_range = [0.15:0.005:0.5]*1e-9;    % list of applied currents

initial_rate = zeros(size(current_range)); % array to store 1/(first ISI)
final_rate = zeros(size(current_range));   % array to store 1/(final ISI)
single_spike = zeros(size(current_range)); % array to store "1" for only 1 spike
mean_V = zeros(size(current_range));

for i = 1:length(current_range)           % loop through applied currents
    I_app_b = current_range(i);
    I(:) = I_app_b;                    % update with new applied current for all time points
    
    v = zeros(size(t_vector));       % initialize voltage array
    v(1) = E_L;                     % set value of initial membrane potential
    I_sra = zeros(size(t_vector));   % initialize adaptation current
    spikes = zeros(size(t_vector));  % initialize vector to store spikes
    
    for j = 1:length(t_vector)-1     % simulation for all time points
        
        if ( v(j) > V_max )              % if there is a spike
            v(j) = V_reset;             % reset the voltage
            I_sra(j) = I_sra(j) + b;    % increase the adaptation current by b
            spikes(j) = 1;              % record the spike
        end
        
        % next line integrates the voltage over time, first part is like LIF
        % second part is an exponential spiking term
        % third part includes adaptation
        v(j+1) = v(j) + dt*( g_L*(E_L-v(j) + Delta_th*exp((v(j)-V_th)/Delta_th) ) ...
            - I_sra(j) + I(j))/C_m;
        
        % next line decays the adaptation toward a steady state in between spikes
        I_sra(j+1) = I_sra(j) + dt*( a*(v(j)-E_L) - I_sra(j) )/tau_SRA;
        
    end
    
    spike_times = dt*find(spikes);           % extract the spike times
    
    if ( length(spike_times) > 1 )           % if there is more than 1 spike
        ISIs = diff(spike_times);            % ISI = interval between spikes
        initial_rate(i) = 1/ISIs(1);         % inverse of first ISI
        if ( length(ISIs) > 1 )             % if there are further ISIs
            final_rate(i) = 1/ISIs(end);     % inverse of final ISI
        end
    else
        if ( length(spike_times) == 1 )      % if there is only one spike
            single_spike(i) = 1;             % record "1" for this trial
        end
    end
    
    mean_V(i) = mean(v);
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