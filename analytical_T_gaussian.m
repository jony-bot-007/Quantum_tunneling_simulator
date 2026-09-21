function T = analytical_T_gaussian(E, V0, x1, width, m, hbar)
T = zeros(size(E));

for i = 1:numel(E)
    Ei = E(i);
    sigma = width / (2 * sqrt(2 * log(2)));
    epsilon = (pi * sigma * sqrt(m)) / hbar * (V0 - Ei) ./ sqrt(V0);
    T = 1 ./ (1 + exp(epsilon));
end
end
