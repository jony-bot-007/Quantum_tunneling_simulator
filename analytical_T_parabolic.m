function T = analytical_T_parabolic(E, V0, x1, width, m, hbar)
T = zeros(size(E));
for i = 1:numel(E)
    Ei = E(i);
    if abs(Ei - V0) < 1e-9
        T(i) = 0.5;
    elseif Ei < V0
        T(i) = 1 ./ (1 + exp( (pi * sqrt(m * width^2) / (hbar * sqrt(2 * V0))) * (V0 - E) ));
    else
        T(i) = 1 ./ (1 + exp( -(pi * sqrt(m * width^2) / (hbar * sqrt(2 * V0))) * (E - V0) ));
    end
end
end
