clear; close all; clc;

% Input parameters
a   = 1e-5;     % m        % Characteristic pore/channel radius
k   = 5e-14;    % m/Pa.s   % Hydraulic permeability coefficient of the medium
T   = 293.15;   % K        % Absolute temperature of the system
Rg  = 8.3145;   % J/mol.K  % Universal gas constant 
l   = 10;       % m        % Total length of the domain
c0  = 800;      % mol/m^3   % Inlet boundary concentration
psil_L= -1e6;    % Pa     % Leaf xylem water potential
psil_s= -0.1e6;
mu  = 1.5e-3;   % Pa.s     % Dynamic viscosity of the fluid
D   = 4e-10;    % m^2/s    % Molecular diffusion coefficient of the species
tr  = 1000;     % s
g = 9.81;    % gravity
rho = 1000;  % density of water

% nondimensional numbers
X=abs(psil_L)/(Rg*T*c0);
M=(16*k*mu*l^2)/a^3;   % Münch number
p_0 = Rg*T*c0; % pressure scale
u_0=(a^2*p_0)/(8*mu*l);
Pe=(u_0*l)/D;          % Peclet number
t_0=l/u_0; % advective timescale
Da=t_0/tr;         %  Damköhler number
ph=rho*g*l; % hydrostatic pressure

% Domain (input or predefined)
n   = 100;
dz  = l / n;
duration = 100;
% dt_CFL = dz^2 / (2*D);% maximum stable time when you have diffusion
% dt = min(100.0, dt_CFL); % to avoid large jumps in  time
dt=0.01;
% xylem water potential
psil = (psil_L - psil_s) .* (1 - (dz/l) .* ((1:n)' - 0.5)) + psil_s; % xylem water potential
Psi=psil./abs(psil_L);

% INITIAL GUESS
c_old_init = (c0 - (rho*g*l - psil(n))/(Rg*T) ).* (1 - (dz/l) .* ((1:n)' - 0.5)) ...
    + (rho*g*l - psil(n))/(Rg*T);   % uniform
C_old_scaled=c_old_init./c0;
dz = dz./l; % scaled
P_old = solve_pressure(C_old_scaled,Psi,X,n,M);


% for saving
C= zeros(n,duration);
P = zeros(n,duration);

C(:,1) = C_old_scaled;
P(:,1) = P_old;


for i = 2:duration

    [P(:,i),C(:,i),~] = Newton(n,P(:,i-1),M,X,Pe,Da,ph,p_0,Psi,C(:,i-1),dt);
    
end

figure
plot(C(:,1))
hold on
plot(C(:,end))


figure
plot(P(:,1))
hold on
plot(P(:,end))