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
tr  = 10000000;    % [s]      Homogeneous reaction time scale
trs = tr/a;    % [s/m]    Surface reaction time scale
g   = 9.81;      % [m/s^2]  Gravity
rho = 1000;      % [kg/m^3] Density of water
J_assim = 1e-6;% [mol m^-2 s^-1] sucrose assimilation flux;


% 2) NON-DIMENSIONAL NUMBERS

X   = abs(psil_L)/(Rg*T*c0);        % Osmotic number
M   = (16*k*mu*l^2)/a^3;            % Münch number
p_0 = Rg*T*c0;                      % Pressure scale
u_0 = (a^2*p_0)/(8*mu*l);           % Characteristic velocity
Pe  = (u_0*l)/D;                    % Peclet number (advection/diffusion)
t_0 = l/u_0;                        % Advective time scale
Da  = t_0/tr;                       % Homogeneous Damköhler number
Das = t_0/(trs*l);                  % Heterogeneous Damköhler number
ph  = (rho*g*l)/p_0;                      % Hydrostatic pressure scale
J_assim = J_assim/(u_0*c0);         %  sucrose assimilation flux;

% 3) NUMERICAL DOMAIN & DISCRETIZATION

n   = 100;          % Number of spatial grid points
dz  = l / n;        % Spatial step (dimensional)
dZ  = dz / l;       % Non-dimensional spatial step
duration = 1000;     % Number of time steps
dt = 0.001;         % Time step size (non-dimensional)
threshold = 0;      % Switch, 0: Cs = 0
                    %         1: Cs = - d2P/dZ2 + ph/p0 - X*Psi 

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

% c0 - (c0 - rho*g*l/(Rg*T) + psil(n)/(Rg*T)).*...
%     (dz/l) .* ((1:n)' - 0.5);
% Non-dimensional concentration
% C_old_scaled = c_old_init ./ c0;
C_old_scaled = ph - 3*X.*Psi;

% Initial pressure from steady-state solver
P_old = solve_pressure(C_old_scaled, Psi, X, n, M, dZ);



% 6) STORAGE MATRICES

C = zeros(n, duration);   % Concentration matrix
P = zeros(n, duration);   % Pressure matrix

C(:,1) = C_old_scaled;    % Store initial concentration
P(:,1) = P_old;           % Store initial pressure;

% 7) TIME INTEGRATION LOOP (Newton Solver)

for i = 2:duration
    
    % Solve nonlinear coupled system using Newton's method
    [P(:,i), C(:,i), ~] =  Newton(n, P(:,i-1), M, X, Pe, Da, Das, ...
               ph, Psi, C(:,i-1), dt, dZ,J_assim,threshold);
    
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




% UPDATED CONCENTRATION GRAPH
% INITIAL + 2 VALUES OF J_assim


figure

% INITIAL PROFILE


plot(C(:,1).*c0,'LineWidth',2)

hold on


% J_assim = 10^-3


J_temp = 1e-6;

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X, Pe, Da, Das, ...
        ph, Psi, C_temp, dt, dZ, J_temp, threshold);

end

plot(C_temp.*c0,'LineWidth',2)


% J_assim = 10^-1


J_temp = 1e-1;

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X, Pe, Da, Das, ...
        ph, Psi, C_temp, dt, dZ, J_temp, threshold);

end

plot(C_temp.*c0,'LineWidth',2)


% GRAPH SETTINGS


title('Effect of Assimilation Flux on concentration')

xlabel('Distance')

ylabel('Concentration [mol/m^3]')

legend('Initial', ...
       'J_{assim}=10^{-6}', ...
       'J_{assim}=10^{-1}')


box on

set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

set(gcf,'Color','w')

exportgraphics(gcf,'assimilation_concentration.png','Resolution',300)


% UPDATED PRESSURE GRAPH
% INITIAL + 2 VALUES OF PSI


figure


% INITIAL PROFILE


plot(P(:,1).*p_0,'LineWidth',2)

hold on


% PSI = -0.2 MPa


psil_L_temp = -0.2e6;

X_temp = abs(psil_L_temp)/(Rg*T*c0);

psil_temp = (psil_L_temp - psil_s) .* ...
    (1 - (dz/l) .* ((1:n)' - 0.5)) + psil_s;

Psi_temp = psil_temp ./ abs(psil_L_temp);

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X_temp, Pe, Da, Das, ...
        ph, Psi_temp, C_temp, dt, dZ, J_assim, threshold);

end

plot(P_temp.*p_0,'LineWidth',2)


% PSI = -1 MPa


psil_L_temp = -1e6;

X_temp = abs(psil_L_temp)/(Rg*T*c0);

psil_temp = (psil_L_temp - psil_s) .* ...
    (1 - (dz/l) .* ((1:n)' - 0.5)) + psil_s;

Psi_temp = psil_temp ./ abs(psil_L_temp);

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X_temp, Pe, Da, Das, ...
        ph, Psi_temp, C_temp, dt, dZ, J_assim, threshold);

end

plot(P_temp.*p_0,'LineWidth',2)


% GRAPH SETTINGS


title('Effect of Xylem Water Potential on pressure')

xlabel('Distance')

ylabel('Pressure [Pa]')

legend('Initial', ...
       '\psi_L=-0.2 MPa', ...
       '\psi_L=-1 MPa')


box on

set(gca,'FontSize',14)
set(gca,'LineWidth',1.5)

set(gcf,'Color','w')



% BETTER PRESSURE GRAPH APPEARANCE


title('Effect of Xylem Water Potential on Pressure')

ylabel('Pressure [Pa]')

legend('Initial', ...
       '\psi_L=-0.2 MPa', ...
       '\psi_L=-1 MPa')
exportgraphics(gcf,'psi_pressure.png','Resolution',300)


% FIX THE SCALE LIKE THE ORIGINAL GRAPH


xlim([1 n])

ylim([0.2e6 2e6])




box on

set(gca,'FontSize',14)

set(gca,'LineWidth',1.5)

set(gcf,'Color','w')


% EXPORT CONCENTRATION GRAPH


exportgraphics(gcf,'assimilation_effect.png','Resolution',300)


% EFFECT OF ASSIMILATION FLUX ON PRESSURE


figure


% INITIAL PRESSURE


plot(P(:,1).*p_0,'LineWidth',2)

hold on


% J_assim = 10^-6


J_temp = 1e-6;

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X, Pe, Da, Das, ...
        ph, Psi, C_temp, dt, dZ, J_temp, threshold);

end

plot(P_temp.*p_0,'LineWidth',2)

% J_assim = 10^-1


J_temp = 1e-1;

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X, Pe, Da, Das, ...
        ph, Psi, C_temp, dt, dZ, J_temp, threshold);

end

plot(P_temp.*p_0,'LineWidth',2)


% GRAPH SETTINGS

title('Effect of Assimilation Flux on Pressure',...
'FontSize',16,...
'FontWeight','bold')

xlabel('Distance')

ylabel('Pressure [Pa]')

legend('Initial', ...
       'J_{assim}=10^{-6}', ...
       'J_{assim}=10^{-1}')

xlim([1 n])

ylim([0.2e6 2e6])

box on

set(gca,'FontSize',14)

set(gca,'LineWidth',1.5)

set(gcf,'Color','w')

exportgraphics(gcf,'assimilation_pressure.png','Resolution',300)

% EFFECT OF XYLEM WATER POTENTIAL ON CONCENTRATION


figure


% INITIAL CONCENTRATION

plot(C(:,1).*c0,'LineWidth',2)

hold on


% PSI = -0.2 MPa


psil_L_temp = -0.2e6;

X_temp = abs(psil_L_temp)/(Rg*T*c0);

psil_temp = (psil_L_temp - psil_s) .* ...
    (1 - (dz/l) .* ((1:n)' - 0.5)) + psil_s;

Psi_temp = psil_temp ./ abs(psil_L_temp);

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X_temp, Pe, Da, Das, ...
        ph, Psi_temp, C_temp, dt, dZ, J_assim, threshold);

end

plot(C_temp.*c0,'LineWidth',2)


% PSI = -1 MPa


psil_L_temp = -1e6;

X_temp = abs(psil_L_temp)/(Rg*T*c0);

psil_temp = (psil_L_temp - psil_s) .* ...
    (1 - (dz/l) .* ((1:n)' - 0.5)) + psil_s;

Psi_temp = psil_temp ./ abs(psil_L_temp);

C_temp = C_old_scaled;
P_temp = P_old;

for i = 2:100

    [P_temp, C_temp, ~] = Newton( ...
        n, P_temp, M, X_temp, Pe, Da, Das, ...
        ph, Psi_temp, C_temp, dt, dZ, J_assim, threshold);

end

plot(C_temp.*c0,'LineWidth',2)


% GRAPH SETTINGS


title('Effect of Xylem Water Potential on Concentration')

xlabel('Distance')

ylabel('Concentration [mol/m^3]')

legend('Initial', ...
       '\psi_L=-0.2 MPa', ...
       '\psi_L=-1 MPa')

xlim([1 n])

ylim([0 1400])

box on

set(gca,'FontSize',14)

set(gca,'LineWidth',1.5)

set(gcf,'Color','w')

exportgraphics(gcf,'psi_concentration.png','Resolution',300)