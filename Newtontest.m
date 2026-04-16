function [J_an, J_num, maxdiff, maxrel] = Newtontest(n,P,C,M,X,Pe,Da,Das,ph,Psi,C_old_scaled,dt,dZ,J_assim,threshold)
%NEWTONTEST Compare analytic Jacobian to finite-difference Jacobian.
%   [J_AN, J_NUM, MAXDIFF, MAXREL] = NEWTONTEST(...) computes the
%   analytic Jacobian from Jacobian(...) and a numerical Jacobian via
%   finite differences on Fvector(...). It returns the two Jacobians and
%   the maximum absolute and relative differences.
%
%   Call it with the same arguments you use for Newton: the current
%   pressure and concentration profiles plus all model parameters.

    % Accept either vectors or storage matrices for P and C.
    if size(P,2) > 1
        P = P(:,1);
    end
    if size(C,2) > 1
        C = C(:,1);
    end

    U = zeros(2*n,1);
    U(1:2:2*n) = P;
    U(2:2:2*n) = C;

    % Compute base residual and analytic Jacobian.
    F0 = Fvector(n,P,C,M,X,Pe,Da,Das,ph,Psi,C_old_scaled,dt,dZ,J_assim,threshold);
    J_an = Jacobian(n,P,C,M,Da,Das,Pe,dt,dZ,threshold);

    % Finite-difference Jacobian.
    J_num = zeros(2*n,2*n);
    eps_rel = 1e-6;

    for k = 1:2*n
        h = eps_rel * max(1, abs(U(k)));
        if h == 0
            h = eps_rel;
        end

        Upert = U;
        Upert(k) = Upert(k) + h;

        Pp = Upert(1:2:2*n);
        Cp = Upert(2:2:2*n);

        Fp = Fvector(n,Pp,Cp,M,X,Pe,Da,Das,ph,Psi,C_old_scaled,dt,dZ,J_assim,threshold);
        J_num(:,k) = (Fp - F0) / h;
    end

    diff = J_an - J_num;
    maxdiff = max(abs(diff), [], 'all');
    maxrel = max(abs(diff) ./ max(1e-12, abs(J_num)), [], 'all');

    fprintf('Newtontest: max abs diff = %g, max rel diff = %g\n', maxdiff, maxrel);
end