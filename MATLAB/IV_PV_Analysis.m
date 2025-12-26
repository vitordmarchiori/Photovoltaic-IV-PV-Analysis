%% ========================================================================
% PV CELL SIMULATOR
% MATLAB IMPLEMENTATION - I-V AND P-V CURVES UNDER VARIOUS CONDITIONS
%
% Vítor D. Marchiori, Rodrigo B. Santos and Douglas D. Bueno
% ========================================================================
clear; close all; clc;

%% ---------------------- PHYSICAL CONSTANTS ------------------------------
k  = 1.380649e-23;       % Boltzmann constant [J/K]
q  = 1.602176634e-19;    % Electron charge [C]

%% ---------------------- CELL / MODULE PARAMETERS ------------------------
Isc_ref    = 5.1;        % Short-circuit current at reference conditions [A]
Voc_ref    = 42;         % Open-circuit voltage at reference conditions [V]
Rs_default = 4e-3;       % Series resistance [Ohm]
Rsh_default= 1000;       % Shunt resistance [Ohm]
Ns         = 36;         % Number of cells in series per module

%% ---------------------- MODEL PARAMETERS --------------------------------
n      = 1.2;            % Diode ideality factor
G_ref  = 1000;           % Reference irradiance [W/m²]
T_ref  = 25;             % Reference temperature [°C]
Ki     = 0.0017;         % Temperature coefficient of Isc [A/°C]
Eg     = 1.12;           % Bandgap energy [eV]
Tnom   = 298;            % Nominal temperature [K] (25°C)
Is_ref = 1e-9;           % Reverse saturation current at Tnom [A]

%% ---------------------- INPUT RANGES FOR PARAMETRIC STUDIES ------------
temperatures = [0 25 50 75 100];            % [°C] – temperature sweep
irradiances  = [200 400 600 800 1000];      % [W/m²] – irradiance sweep
Rs_values    = [1e-3 4e-3 8e-3 12e-3 16e-3];% [Ohm] – series resistance sweep
Rsh_values   = [10 100 400 700 1000];       % [Ohm] – shunt resistance sweep
Is_values    = [1e-6 1e-7 1e-8 1e-9 1e-10]; % [A] – reverse saturation current sweep

%% ----------------------- MODEL FUNCTIONS --------------------------------
% Photocurrent (Iph) as a function of irradiance G [W/m²] and temperature T_C [°C]
% (Equivalent form of the expression used by Selmi et al.)
photo_current = @(G, T_C) (Isc_ref + Ki*(T_C - 298)) .* (G / G_ref);

% Reverse saturation current Is(T) as a function of cell temperature T_K [K]
% Is(T) = Is(Tnom) * (T/Tnom)^3 * exp( (q*Eg/(n*k)) * (1/Tnom - 1/T) )
reverse_saturation_current = @(T_K) ...
    Is_ref .* (T_K./Tnom).^3 .* exp( (q*Eg/(n*k)) .* (1./Tnom - 1./T_K) );

% Thermal voltage per cell including ideality factor:
% Vt_cell = n * k * T / q
Vt_cell = @(T_K) (n*k*T_K/q); 

% Robust I(V) solver (continuation + fallback) with exponential cap
EXP_CAP = 80;  % limits the exponential argument to avoid overflow

solve_I = @(V, Iph, Is, T_K, Rs, Rsh, Iguess) local_solve_I( ...
    V, Iph, Is, T_K, Rs, Rsh, Iguess, Ns, n, k, q, EXP_CAP);

%% ---------------------- VOLTAGE SWEEP -----------------------------------
% Voltage range for I-V and P-V curves
V = linspace(0, 42, 220);

% Helper to standardize axis style (LaTeX labels, grid, etc.)
style_axes = @() set(gca, ...
    'TickLabelInterpreter','latex', ...
    'FontSize',14, ...
    'XMinorGrid','on', ...
    'YMinorGrid','on', ...
    'GridAlpha',0.3);

%% ========================================================================
% 1 - Effect of Irradiance
% ========================================================================
figure('Color','w','Position',[100,100,800,500],'Name','Effect of Irradiance');

% I-V characteristics
subplot(1,2,1); hold on;
xlabel('Voltage [V]','Interpreter','latex','FontSize',14);
ylabel('Current [A]','Interpreter','latex','FontSize',14);
grid on; box on; style_axes();
ylim([0 6]); xlim([0 26]);

% P-V characteristics
subplot(1,2,2); hold on;
xlabel('Voltage [V]','Interpreter','latex','FontSize',14);
ylabel('Power [W]','Interpreter','latex','FontSize',14);
grid on; box on; style_axes();
ylim([0 100]); xlim([0 26]);

for G = irradiances
    T_C = T_ref;
    T_K = T_C + 273.15;
    Iph = photo_current(G, T_C);
    Is  = reverse_saturation_current(T_K);

    I = zeros(size(V));
    Iguess = Iph;  % initial guess close to short-circuit current

    for j = 1:numel(V)
        I(j) = solve_I(V(j), Iph, Is, T_K, Rs_default, Rsh_default, Iguess);
        % continuation: use the previous solution as the next initial guess
        Iguess = max(0, min(Iph, I(j)));
    end

    P = V .* I;

    subplot(1,2,1);
    plot(V, I, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d W/m$^2$', G));

    subplot(1,2,2);
    plot(V, P, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d W/m$^2$', G));
end

subplot(1,2,1);
legend('show','Interpreter','latex','Location','northeast');

subplot(1,2,2);
legend('show','Interpreter','latex','Location','northwest');


%% ========================================================================
% 2 - Effect of Temperature
% ========================================================================
figure('Color','w','Position',[100,100,800,500],'Name','Effect of Temperature');

% I-V characteristics
subplot(1,2,1); hold on;
xlabel('Voltage [V]','Interpreter','latex','FontSize',14);
ylabel('Current [A]','Interpreter','latex','FontSize',14);
grid on; box on; style_axes();
ylim([0 5]); xlim([0 28]);

% P-V characteristics
subplot(1,2,2); hold on;
xlabel('Voltage [V]','Interpreter','latex','FontSize',14);
ylabel('Power [W]','Interpreter','latex','FontSize',14);
grid on; box on; style_axes();
ylim([0 105]); xlim([0 28]);

for T_C = temperatures
    T_K = T_C + 273.15;
    Iph = photo_current(G_ref, T_C);
    Is  = reverse_saturation_current(T_K);

    I = zeros(size(V));
    Iguess = Iph;

    for j = 1:numel(V)
        I(j) = solve_I(V(j), Iph, Is, T_K, Rs_default, Rsh_default, Iguess);
        Iguess = max(0, min(Iph, I(j)));
    end

    P = V .* I;

    subplot(1,2,1);
    plot(V, I, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d$^\\circ$C', T_C));

    subplot(1,2,2);
    plot(V, P, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d$^\\circ$C', T_C));
end

subplot(1,2,1);
legend('show','Interpreter','latex','Location','southwest');

subplot(1,2,2);
legend('show','Interpreter','latex','Location','northwest');

%% ========================================================================
% 3 - Effect of Series Resistance Rs
% ========================================================================
figure('Color','w','Position',[100,100,800,500],'Name','Effect of Rs');

% I-V characteristics
subplot(1,2,1); hold on;
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Current [A]','Interpreter','latex');
grid on; box on; style_axes();
ylim([0 5]); xlim([0 26]);

% P-V characteristics
subplot(1,2,2); hold on;
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Power [W]','Interpreter','latex');
grid on; box on; style_axes();
ylim([0 96]); xlim([0 26]);

T_C = T_ref;
T_K = T_C + 273.15;
Iph = photo_current(G_ref, T_C);
Is  = reverse_saturation_current(T_K);

for Rs = Rs_values
    I = zeros(size(V));
    Iguess = Iph;

    for j = 1:numel(V)
        I(j) = solve_I(V(j), Iph, Is, T_K, Rs, Rsh_default, Iguess);
        Iguess = max(0, min(Iph, I(j)));
    end

    P = V .* I;

    % Label using standard scientific notation
    label_str = sprintf('%.1e $\\Omega$', Rs);

    subplot(1,2,1);
    plot(V, I, 'LineWidth', 2, 'DisplayName', label_str);

    subplot(1,2,2);
    plot(V, P, 'LineWidth', 2, 'DisplayName', label_str);
end

subplot(1,2,1);
legend('show','Interpreter','latex','Location','southwest');

subplot(1,2,2);
legend('show','Interpreter','latex','Location','northwest');

%% ========================================================================
% 4 - Effect of Shunt Resistance Rsh
% ========================================================================
figure('Color','w','Position',[100,100,800,500],'Name','Effect of Rsh');

% I-V characteristics
subplot(1,2,1); hold on;
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Current [A]','Interpreter','latex');
grid on; box on; style_axes();
ylim([0 5]); xlim([0 26]);

% P-V characteristics
subplot(1,2,2); hold on;
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Power [W]','Interpreter','latex');
grid on; box on; style_axes();
ylim([0 96]); xlim([0 26]);

T_C = T_ref;
T_K = T_C + 273.15;
Iph = photo_current(G_ref, T_C);
Is  = reverse_saturation_current(T_K);

for Rsh = Rsh_values
    I = zeros(size(V));
    Iguess = Iph;

    for j = 1:numel(V)
        I(j) = solve_I(V(j), Iph, Is, T_K, Rs_default, Rsh, Iguess);
        Iguess = max(0, min(Iph, I(j)));
    end

    P = V .* I;

    subplot(1,2,1);
    plot(V, I, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d $\\Omega$', Rsh));

    subplot(1,2,2);
    plot(V, P, 'LineWidth', 2, ...
        'DisplayName', sprintf('%d $\\Omega$', Rsh));
end

subplot(1,2,1);
legend('show','Interpreter','latex','Location','southwest');

subplot(1,2,2);
legend('show','Interpreter','latex','Location','northwest');

%% ========================================================================
% 5 - Effect of Reverse Saturation Current Is
% ========================================================================
figure('Color','w','Position',[100,100,800,500],'Name','Effect of Is');

% I-V characteristics
subplot(1,2,1); hold on;
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Current [A]','Interpreter','latex');
grid on; box on; style_axes();
ylim([0 5]); xlim([0 30]);

% P-V characteristics
subplot(1,2,2); hold on;
xlabel('Voltage [V]','Interpreter','latex');
ylabel('Power [W]','Interpreter','latex');
grid on; box on; style_axes();
ylim([0 110]); xlim([0 30]);

T_C = T_ref;
T_K = T_C + 273.15;
Iph = photo_current(G_ref, T_C);

for Is_test = Is_values
    I = zeros(size(V));
    Iguess = Iph;

    for j = 1:numel(V)
        I(j) = solve_I(V(j), Iph, Is_test, T_K, Rs_default, Rsh_default, Iguess);
        Iguess = max(0, min(Iph, I(j)));
    end

    P = V .* I;

    % Legend label in standard scientific notation
    label_str = sprintf('%.1e A', Is_test);

    subplot(1,2,1);
    plot(V, I, 'LineWidth', 2, 'DisplayName', label_str);

    subplot(1,2,2);
    plot(V, P, 'LineWidth', 2, 'DisplayName', label_str);
end

subplot(1,2,1);
legend('show','Interpreter','latex','Location','southwest');

subplot(1,2,2);
legend('show','Interpreter','latex','Location','northwest');

fprintf('All simulations completed successfully.\n');

%% ========================================================================
% LOCAL FUNCTION: NUMERICAL SOLVER FOR THE I–V EQUATION
% ========================================================================
function I = local_solve_I(V, Iph, Is, T_K, Rs, Rsh, Iguess, Ns, n, k, q, EXP_CAP)
%LOCAL_SOLVE_I  Solve the single-diode model equation for current I(V).
%
%   Given:
%       V      - terminal voltage [V]
%       Iph    - photocurrent [A]
%       Is     - reverse saturation current [A]
%       T_K    - cell temperature [K]
%       Rs     - series resistance [Ohm]
%       Rsh    - shunt resistance [Ohm]
%       Iguess - initial guess for current [A]
%       Ns     - number of series cells
%       n      - diode ideality factor
%       k, q   - physical constants
%       EXP_CAP- upper limit for exponential argument to avoid overflow
%
%   The function uses fzero with a continuation strategy and basic
%   bracketing logic to improve robustness.

    % Nonlinear I-V equation in implicit form f(I) = 0
    f = @(I) I - Iph ...
        + Is .* ( exp( min( q*(V + I.*Rs) / (Ns*n*k*T_K), EXP_CAP) ) - 1 ) ...
        + (V + I.*Rs) ./ Rsh;

    % First attempt: single fzero call using the given guess
    try
        I = fzero(f, Iguess);
        if isfinite(I), return; end
    catch
        % If fzero fails or throws, we continue and try to bracket
    end

    % Bracketing strategy in case the first attempt failed:
    % start from a lower bound close to zero and an upper bound near Iph
    Il = max(0, min(Iguess, 0.2*Iph));
    Iu = max(Iph, Iguess + 1e-6);

    fl = f(Il);
    fu = f(Iu);

    it = 0;
    ok = isfinite(fl) && isfinite(fu) && sign(fl) * sign(fu) <= 0;

    % Expand the bracket until a valid sign change is found or max iterations reached
    while ~ok && it < 20
        Il = max(0, 0.5 * Il);
        Iu = Iu * 1.5 + 1e-6;

        fl = f(Il);
        fu = f(Iu);

        ok = isfinite(fl) && isfinite(fu) && sign(fl) * sign(fu) <= 0;
        it = it + 1;
    end

    % If a valid bracket exists, use fzero on the interval [Il, Iu]
    if ok
        I = fzero(f, [Il, Iu]);
    else
        % Fallback: clamp the value within [0, Iph] using the last guess
        I = max(0, min(Iph, Iguess));
    end
end
