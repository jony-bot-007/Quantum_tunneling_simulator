function T = analytical_T_square(E, V0, L, m, hbar)
T = zeros(size(E));
for i = 1:numel(E)
    Ei = E(i);
    if abs(Ei - V0) < 1e-9
        T(i) = 1 / (1 + m*V0*L^2 / (2*hbar^2));
    elseif Ei < V0
        kappa = sqrt(2*m*(V0 - Ei)) / hbar;
        T(i)  = 1 / (1 + (V0^2 * sinh(kappa*L)^2) / (4*Ei*(V0 - Ei)));
    else
        k2   = sqrt(2*m*(Ei - V0)) / hbar;
        T(i) = 1 / (1 + (V0^2 * sin(k2*L)^2) / (4*Ei*(Ei - V0)));
    end
end
end
