function T = analytical_T_trapezoidal(E, V0, x1, width, m, hbar)
T = zeros(size(E));
b = 0.50*width;
d = 0.25*width;
for i = 1:numel(E)
    Ei = E(i);
     if abs(Ei - V0) < 1e-9
        T = 1 / (1 + (m * V0 * (b + 0.67 * d)^2) / (2 * hbar^2));
    elseif Ei < V0
        kappa = sqrt(2 * m * (V0 - E)) / hbar;
        gamma = (4 * sqrt(2 * m) ./ (3 * hbar * V0)) * d .* ((V0 - E).^(1.5)) + 2 * kappa * b;
        T = exp(-gamma);

     else
        k = sqrt(2 * m * E) / hbar;
        k1 = sqrt(2 * m * (E - V0)) / hbar;
        T = (4 * k .* k1) ./ ((k + k1).^2);
     end
end
end
