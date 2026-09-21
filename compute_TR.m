function [T, R, prob_total] = compute_TR(x, psi, barrier_x0, barrier_x1)
dx = x(2) - x(1);
prob_density = abs(psi).^2;
prob_total = sum(prob_density) * dx;
T = sum(prob_density(x > barrier_x1)) * dx / prob_total;
R = sum(prob_density(x < barrier_x0)) * dx / prob_total;
end
