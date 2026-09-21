function T = analytical_T_double_square(E, V0, x1, width, m, hbar)
T = zeros(size(E));

a = 0.25 * width;      % width of each barrier
L = 0.50 * width;      % gap between barriers

for i = 1:numel(E)
    Ei = E(i);
    if Ei < V0
        k = sqrt(2*m*Ei) / hbar;
        kappa = sqrt(2*m*(V0-Ei)) / hbar;
        A = cosh(kappa*a) + 1i*(k^2-kappa^2)/(2*k*kappa)*sinh(kappa*a);
        B = 1i*(k^2+kappa^2)/(2*k*kappa)*sinh(kappa*a);
        C = -1i*(k^2+kappa^2)/(2*k*kappa)*sinh(kappa*a);
        D = cosh(kappa*a) - 1i*(k^2-kappa^2)/(2*k*kappa)*sinh(kappa*a);
        M11 = A*exp(1i*k*L)*A + B*exp(-1i*k*L)*C;
        T(i) = 1 / abs(M11)^2;
    elseif abs(Ei-V0) < 1e-9
        k = sqrt(2*m*Ei) / hbar;
        A = 1 + 1i*k*a/2;
        B = 1i*k*a/2;
        C = -1i*k*a/2;
        D = 1 - 1i*k*a/2;
        M11 = A*exp(1i*k*L)*A + B*exp(-1i*k*L)*C;
        T(i) = 1 / abs(M11)^2;
    else
        k = sqrt(2*m*Ei) / hbar;
        k1 = sqrt(2*m*(Ei-V0)) / hbar;
        A = cos(k1*a) + 1i*(k^2+k1^2)/(2*k*k1)*sin(k1*a);
        B = 1i*(k^2-k1^2)/(2*k*k1)*sin(k1*a);
        C = -1i*(k^2-k1^2)/(2*k*k1)*sin(k1*a);
        D = cos(k1*a) - 1i*(k^2+k1^2)/(2*k*k1)*sin(k1*a);
        M11 = A*exp(1i*k*L)*A + B*exp(-1i*k*L)*C;
        T(i) = 1 / abs(M11)^2;
    end
end
end
