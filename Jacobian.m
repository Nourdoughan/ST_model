function J = Jacobian(n,P,C,M,Da,Das,Pe,dt,dZ)


% JACOBIAN MATRIX FOR COUPLED PRESSURE–CONCENTRATION SYSTEM

% Constructs the 2n × 2n Jacobian matrix associated with the
% nonlinear residual F(P,C).

% Unknown ordering:
% U = [P1,C1,P2,C2,...,Pn,Cn]^T

% Odd indices  -> Pressure equations
% Even indices -> Concentration equations


% n     : number of spatial nodes
% P     : current pressure vector (non-dimensional)
% C     : current concentration vector (non-dimensional)
% M     : Münch number
% Da    : homogeneous Damköhler number
% Das   : heterogeneous Damköhler number
% Pe    : Peclet number
% dt    : time step
% dZ    : non-dimensional spatial step

% Internal node indices
i = 2:n-1; % location vector

% Allocate full Jacobian matrix
J = zeros(2*n,2*n);


% PRESSURE EQUATION (internal nodes)
% Residual form:
% (-1/(M dZ^2))*(P(i+1)-2P(i)+P(i-1)) + P(i) - C(i)

% ∂Fp/∂P(i-1)
J(sub2ind(size(J),2.*i-1,2.*(i-1)-1)) = -1/(M*dZ^2);

% ∂Fp/∂P(i)
J(sub2ind(size(J),2.*i-1,2.*i-1)) = 2/(M*dZ^2) + 1;

% ∂Fp/∂C(i)
J(sub2ind(size(J),2.*i-1,2.*i)) = -1;

% ∂Fp/∂P(i+1)
J(sub2ind(size(J),2.*i-1,2.*(i+1)-1)) = -1/(M*dZ^2);


% CONCENTRATION EQUATION (internal nodes) 
% Residual includes:
% - time derivative
% - advection (pressure-driven)
% - diffusion
% - homogeneous reaction

% ∂Fc/∂P(i-1)
J(sub2ind(size(J),2.*i,2.*(i-1)-1)) = (1/dZ^2).*C(1:n-2);

% ∂Fc/∂C(i-1)
J(sub2ind(size(J),2.*i,2.*(i-1))) = ...
    1/dZ^2.*( P(2:n-1)-P(1:n-2) ) - 1/(Pe*dZ^2);

% ∂Fc/∂P(i)
J(sub2ind(size(J),2.*i,2.*i-1)) = ...
    (1/dZ^2).*( C(2:n-1) + C(1:n-2) );

% ∂Fc/∂C(i)
J(sub2ind(size(J),2.*i,2.*i)) = ...
    1/dt ...
    - 1/dZ^2.*( P(3:n) - P(2:n-1) ) ...
    + 2/(Pe*dZ^2) ...
    + Da;

% ∂Fc/∂P(i+1)
J(sub2ind(size(J),2.*i,2.*(i+1)-1)) = (-1/dZ^2).*C(2:n-1);

% ∂Fc/∂C(i+1)
J(sub2ind(size(J),2.*i,2.*(i+1))) = -1/(Pe*dZ^2);


%  BOUNDARY CONDITIONS

% Top boundary (z = 0) 
% Pressure: P(1) - C(1) - X*Psi(1) = 0
J(1,1) = 1;
J(1,2) = -1;

% Concentration: C(1) = 1 (Dirichlet)
J(2,2) = 1;


% Bottom boundary (z = L) 
% Pressure: P(n) - C(n) - X*Psi(n) = 0
J(2*n-1,2*n-1) = 1;
J(2*n-1,2*n)   = -1;

% Concentration:
% Includes time derivative, surface reaction (Das),
% advection, diffusion, and homogeneous reaction

J(2*n,2*n) = ...
    1/dt ...
    + Das/dZ ...
    + (1/dZ^2)*( ( P(n)-P(n-1) ) + (1/Pe) ) ...
    + Da;

% J(2,2)=)=1/dt - 1/dZ^2.*(P(2)-P(1)) +(1/(Pe*dZ^2)).*C(2) + Da;

% J(2*n,2*n)=1;

end