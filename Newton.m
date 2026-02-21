function [P,C,err] =Newton(n,p_old,M,X,Pe,Da,ph,p_0,Psi,C_old_scaled,dt,dZ)
%dz is the grid siz
% a is the phloem thickness in m                          
%k    % m/Pa.s   % Hydraulic permeability  of the membrane
%T      % K        % Absolute temperature of the system
%Rg     % J/mol.K  % Universal gas constant 
%c0    % mol/L    % Inlet boundary concentration
%psil     % Pa              % xylem water potential
%mu    % Pa.s     % Dynamic viscosity of the fluid
% D      % m^2/s   % Molecular diffusion coefficient of the species
%n                    %number of nodes
% dt = min(1.0, dt_CFL);    % to avoid large jumps in  time
% tr is removal rate


% ITERATION SETTINGS
tol = 1e-6;
max_iter = 1000;
err = 1;
iter = 0;
C=C_old_scaled;
P=p_old;

%ITERATIVE LOOP 
while err > tol && iter < max_iter
    iter = iter + 1;

    % save old
    c_IT = C;
    p_IT = P;

    F = Fvector(n,p_IT,c_IT,M,X,Pe,Da,ph,p_0,Psi,C_old_scaled,dt,dZ);

    J = Jacobian(n,p_IT,c_IT,M,Da,Pe,dt,dZ);

    
   
    % solve using tridiagonal block Thomas algorithm
    
    % BLOCK THOMAS 
    sol = blockThomas2x2(J,-F,n);

    % split corrections
    dp = sol(1:2:2*n);
    dc = sol(2:2:2*n);

    % update
    P= p_IT + dp;
    C= c_IT + dc;

    % Check convergence
    err = norm(Fvector(n,p_IT,c_IT,M,X,Pe,Da,ph,p_0,Psi,C_old_scaled,dt,dZ));
    fprintf("Iteration %d   Error = %e\n", iter, err);

end


fprintf("\nCONVERGED after %d iterations. Final error = %.4e\n", iter, err);


end