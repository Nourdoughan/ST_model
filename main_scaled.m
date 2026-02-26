% Simulation – Coupled Pressure & Concentration
% Time-dependent solution using Newton's method
 
clear all; 
close all; 
clc;

% 1) PHYSICAL PARAMETERS (Dimensional)

a   = 1e-5;      % [m]      Characteristic pore/channel radius
k   = 5e-14;     % [m/Pa.s] Hydraulic permeability of the medium
T   = 293.15;    % [K]      Absolute temperature
Rg  = 8.3145;    % [J/mol.K] Universal gas constant
l   = 10;        % [m]      Total length of the domain
c0  = 800;       % [mol/m^3] Inlet boundary concentration
psil_L = -1e6;   % [Pa]     Leaf xylem water potential
psil_s = -0.1e6; % [Pa]     Soil xylem water potential
mu  = 1.5e-3;    % [Pa.s]   Dynamic viscosity
D   = 4e-10;     % [m^2/s]  Molecular diffusion coefficient
tr  = 1000;      % [s]      Homogeneous reaction time scale
trs = 1000/a;    % [s/m]    Surface reaction time scale
g   = 9.81;      % [m/s^2]  Gravity
rho = 1000;      % [kg/m^3] Density of water

% 2) NON-DIMENSIONAL NUMBERS

X   = abs(psil_L)/(Rg*T*c0);        % Osmotic number
M   = (16*k*mu*l^2)/a^3;            % Münch number
p_0 = Rg*T*c0;                      % Pressure scale
u_0 = (a^2*p_0)/(8*mu*l);           % Characteristic velocity
Pe  = (u_0*l)/D;                    % Peclet number (advection/diffusion)
t_0 = l/u_0;                        % Advective time scale
Da  = t_0/tr;                       % Homogeneous Damköhler number
Das = t_0/(trs*l);                  % Heterogeneous Damköhler number
ph  = rho*g*l;                      % Hydrostatic pressure scale

% 3) NUMERICAL DOMAIN & DISCRETIZATION

n   = 100;          % Number of spatial grid points
dz  = l / n;        % Spatial step (dimensional)
dZ  = dz / l;       % Non-dimensional spatial step
duration = 100;     % Number of time steps
dt = 0.001;         % Time step size (non-dimensional)

% 4) XYLEM WATER POTENTIAL (Linear Profile)

% Linear variation between leaf and soil
psil = (psil_L - psil_s) .* ...
       (1 - (dz/l) .* ((1:n)' - 0.5)) + psil_s;

% Non-dimensional xylem potential
Psi = psil ./ abs(psil_L);

% 5) INITIAL CONDITIONS

% Initial concentration profile (hydrostatic correction included)
c_old_init = (c0 - (rho*g*l - psil(n))/(Rg*T)) .* ...
             (1 - (dz/l) .* ((1:n)' - 0.5)) ...
             + (rho*g*l - psil(n))/(Rg*T);

% Non-dimensional concentration
C_old_scaled = c_old_init ./ c0;

% Initial pressure from steady-state solver
P_old = solve_pressure(C_old_scaled, Psi, X, n, M, dZ);



% 6) STORAGE MATRICES

C = zeros(n, duration);   % Concentration matrix
P = zeros(n, duration);   % Pressure matrix

C(:,1) = C_old_scaled;    % Store initial concentration
P(:,1) = P_old;           % Store initial pressure


% 7) TIME INTEGRATION LOOP (Newton Solver)

for i = 2:duration
    
    % Solve nonlinear coupled system using Newton's method
    [P(:,i), C(:,i), ~] =  Newton(n, P(:,i-1), M, X, Pe, Da, Das, ...
               ph, p_0, Psi, C(:,i-1), dt, dZ);
    
end


% 8) POST-PROCESSING & VISUALIZATION

%  Concentration evolution
figure
plot(C(:,1).*c0)           % Initial concentration
hold on
plot(C(:,end).*c0)         % Final concentration
title('Concentration Evolution')
xlabel('Spatial Node')
ylabel('Concentration [mol/m^3]')
legend('Initial','Final')


%  Pressure evolution 
figure
plot(P(:,1).*p_0)          % Initial pressure
hold on
plot(P(:,end).*p_0)        % Final pressure

% Add hydrostatic correction
plot(P(:,end).*p_0 + rho*g.*(dz .* ((1:n)' - 0.5) - l))

title('Pressure Evolution')
xlabel('Spatial Node')
ylabel('Pressure [Pa]')
legend('Initial','Final','Final + Hydrostatic')
