function T = analytical_T_triangular(E, V0, x1, width, m, hbar)
T = zeros(size(E));
for i = 1:numel(E)
    Ei = E(i);
    if abs(Ei - V0) < 1e-9
        T(i) = 1 / (1 + 0.67*((m*V0*width^2) / hbar^2)^(1/3));
    elseif Ei < V0
        T(i)  = exp(-(4*sqrt(2*m)*(V0-E).^(1.5)*width)/(3*hbar*V0));
    else
        k2   = sqrt(2*m*(Ei - V0)) / hbar;
        T(i) = (4*sqrt(E.*(E-V0)))./(sqrt(E)+sqrt(E-V0)).^2;
    end
end
end
