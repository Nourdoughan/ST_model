function F = Fvector(n,P,C,M,X,Pe,Da,ph,p_0,Psi,C_old_scaled,dt,dZ)

%  parameters
% n is the number of nodes
% dz is spatial step size
% a is Pore/channel radius [m]
% mu is Dynamic viscosity of fluid [Pa·s]
% D is Molecular diffusivity of species in fluid [m^2/s]
% k is Hydraulic permeability of porous medium m/(Pa·s) 
% R is Universal gas constant [J/(mol·K)]
% T is Absolute temperature [K]
% psil is xylem water potential
% p = 
% c = 



% allocate F 
F = zeros(2*n,1);

% location vector
i = 2:n-1;

%  p-equation
F(2.*i-1) = (-1/(M*dZ^2)).*( P(i+1) - 2.*P(i) + P(i-1) ) ...
             + P(i) - C(i) - X.*Psi(i);

%  c-equation
F(2.*i) = ( C(i) - C_old_scaled(i) )./dt - 1/dZ^2*(( C(i).*(P(i+1)-P(i)) - C(i-1).*(P(i)-P(i-1))) ) ...
             -(1/(Pe*dZ^2)).*( C(i+1) - 2.*C(i) + C(i-1) ) + Da.*C(i);


% BC at z = 0
F(1) = P(1) - C(1) - X*Psi(1); % water potential equilibirum
F(2) = C(1) - 1; % cts c at the top

% BC at z = L
F(2*n-1) = P(n) - C(n) - X*Psi(n); % water potential equilibirum
F(2*n)   = C(n) - (ph/p_0 - X*Psi(n)); % p = 0 at the bottom

end