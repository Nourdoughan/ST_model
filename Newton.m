function [P,C,err] = Newton(n,p_old,M,X,Pe,Da,Das,ph,p_0,Psi,C_old_scaled,dt,dZ)


% %  NEWTON SOLVER FOR COUPLED PRESSURE–CONCENTRATION SYSTEM

%  Solves the nonlinear system arising from:
%  - Pressure equation (flow driven by osmotic + hydrostatic effects)
%  - Advection–diffusion–reaction equation for concentration
%  The system is solved using Newton–Raphson iterations
%  with a 2x2 block tridiagonal Thomas algorithm



% % INPUT PARAMETERS

% n              : number of spatial nodes
% p_old          : pressure at previous time step (initial guess)
% C_old_scaled   : concentration at previous time step (scaled)

% M              : Münch number
% X              : osmotic number
% Pe             : Peclet number
% Da             : homogeneous Damköhler number
% Das            : heterogeneous Damköhler number

% ph             : hydrostatic pressure scale
% p_0            : pressure scale
% Psi            : non-dimensional xylem water potential

% dt             : time step
% dZ             : non-dimensional spatial step

% % OUTPUT:
% P              : updated pressure profile
% C              : updated concentration profile
% err            : final residual norm

 % 1) ITERATION SETTINGS


tol = 1e-6;        % Convergence tolerance
max_iter = 1000;   % Maximum Newton iterations

err = 1;           % Initial error (start large)
iter = 0;          % Iteration counter

% Initialize solution vectors with previous time step
C = C_old_scaled;
P = p_old;



% 2) NEWTON ITERATION LOOP

while err > tol && iter < max_iter
    
    iter = iter + 1;

    % Store current iterate (used to compute correction)
    c_IT = C;
    p_IT = P;

    
    % 2.1 Compute Residual Vector
   
    % F contains:
    % - Pressure residuals
    % - Concentration residuals
    % Arranged in interleaved form [P1,C1,P2,C2,...]
    F = Fvector(n,p_IT,c_IT,M,X,Pe,Da,Das,...
                ph,p_0,Psi,C_old_scaled,dt,dZ);

    
    % 2.2 Compute Jacobian Matrix
    
    % J is the block tridiagonal Jacobian of the coupled system
    J = Jacobian(n,p_IT,c_IT,M,Da,Das,Pe,dt,dZ);

    
    % 2.3 Solve Linear System for Newton Correction
    
    % J * delta = -F
    % Solved using a specialized 2x2 block Thomas algorithm
    sol = blockThomas2x2(J, -F, n);

    % Extract pressure and concentration corrections
    dp = sol(1:2:2*n);   % Pressure corrections
    dc = sol(2:2:2*n);   % Concentration corrections

    
    % 2.4 Update Solution
    
    P = p_IT + dp;
    C = c_IT + dc;

    
    % 2.5 Convergence Check
 
    % Compute residual norm to monitor convergence
    err = norm(Fvector(n,p_IT,c_IT,M,X,Pe,Da,Das,...
                       ph,p_0,Psi,C_old_scaled,dt,dZ));

    fprintf("Iteration %d   Error = %e\n", iter, err);

end

% 3) FINAL STATUS


fprintf("\nCONVERGED after %d iterations. Final error = %.4e\n", iter, err);

end