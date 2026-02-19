function P = solve_pressure(C_old_scaled,Psi,X,n,M)

    % Build RHS b_p
    b_p = C_old_scaled + X*Psi;

    
   % Build tridiagonal system
    diag_p = ones(n,1) .* (2/M + 1);   %  M is Münch number
    diag_u = ones(n-1,1) .* (-1/M);
    diag_l = diag_u;

    % Boundary conditions
    diag_p(1) = 1; 
    diag_p(n) = 1;
    diag_u(1) = 0;
    diag_l(end) = 0;

    % Solve using Thomas
    P = thomas(diag_l, diag_p, diag_u, b_p);

end