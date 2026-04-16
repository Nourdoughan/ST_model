function F = Fvector(n,P,C,M,X,Pe,Da,Das,ph,Psi,C_old_scaled,dt,dZ,J_assim,threshold)


%  RESIDUAL VECTOR FOR COUPLED PRESSURE–CONCENTRATION SYSTEM
%  Constructs the nonlinear residual F such that:
%        F(U) = 0
%  where U = [P1,C1,P2,C2,...,Pn,Cn]^T
%  Odd entries  -> Pressure equations
%  Even entries -> Concentration equations

% INPUT DESCRIPTION

% n              : number of spatial nodes
% P              : current pressure vector (non-dimensional)
% C              : current concentration vector (non-dimensional)
% M              : Münch number
% X              : osmotic number
% Pe             : Peclet number
% Da             : homogeneous Damköhler number
% Das            : heterogeneous Damköhler number
% ph             : hydrostatic pressure scale
% Psi            : non-dimensional xylem water potential
% C_old_scaled   : concentration at previous time step
% dt             : time step
% dZ             : non-dimensional spatial step
% threshold      : Switch case for Cs (minimum concentration in the phloem)

% 1) ALLOCATE RESIDUAL VECTOR

F = zeros(2*n,1);   % Total unknowns = 2n (P and C at each node)

% 2) INTERNAL NODES (i = 2 ... n-1)

i = 2:n-1;          % location vector

% 2.1 PRESSURE EQUATION (Poisson-type equation)

% - First term  : axial pressure diffusion (Darcy flow)
% - Second term : pressure
% - Third term  : osmotic contribution from concentration
% - Fourth term : imposed xylem water potential


F(2.*i-1) = (-1/(M*dZ^2)).*( P(i+1) - 2.*P(i) + P(i-1) ) ...
            + P(i) - C(i) - X.*Psi(i);



% 2.2 CONCENTRATION EQUATION

% Terms:
% 1) Time derivative (implicit Euler)
% 2) Advection (pressure-driven transport)
% 3) Diffusion
% 4) Homogeneous reaction (Da)

switch threshold

    case 0
        F(2.*i) = ( C(i) - C_old_scaled(i) )./dt ...
          - (1/dZ^2)*( C(i).*(P(i+1)-P(i)) ...
          - C(i-1).*(P(i)-P(i-1)) ) ...
          - (1/(Pe*dZ^2)).*( C(i+1) - 2.*C(i) + C(i-1) ) ...
          + Da.*C(i);
    case 1
        F(2.*i) = ( C(i) - C_old_scaled(i) )./dt ...
          - (1/dZ^2)*( C(i).*(P(i+1)-P(i)) ...
          - C(i-1).*(P(i)-P(i-1)) ) ...
          - (1/(Pe*dZ^2)).*( C(i+1) - 2.*C(i) + C(i-1) ) ...
          + Da.*( C(i) + ( P(i+1) - 2.*P(i) + P(i-1) )./(M*dZ^2) ...
          - ph + X.*Psi(i) ) ;
end




% 3) BOUNDARY CONDITIONS

% 3.1 TOP BOUNDARY (z = 0)


% Axial velocity equal to zero: dP/dZ = 0
F(1) = - ( P(2) - P(1) )/(M* ...
    dZ^2) + P(1) - C(1) - X*Psi(1);

% Concentration: Dirichlet condition (C = 1 at inlet)

switch threshold

    case 0
        F(2) = ( C(1) - C_old_scaled(1) )./dt ...
         - (1/dZ^2)*( C(1).*( P(2)- P(1) ) ...
         + (1/Pe).*( C(2) - C(1) ) ) - J_assim/dZ ...
         + Da.*C(1);
    case 1
        F(2) = ( C(1) - C_old_scaled(1) )./dt ...
         - (1/dZ^2)*( C(1).*( P(2)- P(1) ) ...
         + (1/Pe).*( C(2) - C(1) ) ) -J_assim/dZ ...
         + Da.*( C(1) + ( P(2) - P(1) )/(M*dZ^2) ...
         - ph + X*Psi(1) );
end



% 3.2 BOTTOM BOUNDARY (z = L)

% Pressure: water potential equilibrium
F(2*n-1) = P(n) - C(n) - X*Psi(n);

% Concentration:
% Includes:
% - Time derivative
% - Surface reaction (Das)
% - Advection
% - Diffusion
% - Homogeneous reaction (Da)

switch threshold
    case 0
        F(2*n) = ( C(n) - C_old_scaled(n) )./dt ...
         + ( C(n) + X*Psi(n) - ph )*Das/dZ ...
         + (1/dZ^2)*( C(n-1).*( P(n)-P(n-1) ) ...
         + (1/Pe).*( C(n) - C(n-1) ) ) ...
         + Da.*C(n);
    case 1
        F(2*n) = ( C(n) - C_old_scaled(n) )./dt ...
         + ( C(n) + X*Psi(n) - ph )*Das/dZ ...
         + (1/dZ^2)*( C(n-1).*( P(n)-P(n-1) ) ...
         + (1/Pe).*( C(n) - C(n-1) ) ) ...
         + Da.*( C(n) - ph + X*Psi(n) );
end


end