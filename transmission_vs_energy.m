clc; clear all;

p.hbar = 1;
p.m    = 1;

p.xmin = -80;
p.xmax = 80;
p.N    = 1000;

p.x0    = -20;
p.sigma = 6;

p.V0            = 1.5;
p.barrier_width = 4;
p.barrier_x0    = -p.barrier_width/2;

p.barrier_type  = 'double_square'; % Any of 'square', 'triangular', 'gaussian', 'parabolic', 'trapezoidal', 'double_square', 'double_gaussian'.

p.dt            = 0.02;
p.nsteps        = 3000;
p.store_history = false;

x = linspace(p.xmin, p.xmax, p.N)';
V = build_potential(x, p);

switch lower(p.barrier_type)

    case 'square'
        T_of_E = @(Ej) analytical_T_square(Ej, p.V0, p.barrier_width, p.m, p.hbar);

    case 'triangular'
        T_of_E = @(Ej) analytical_T_triangular(Ej, p.V0, p.barrier_x0, p.barrier_width, p.m, p.hbar);

    case 'parabolic'
        T_of_E = @(Ej) analytical_T_parabolic(Ej, p.V0, p.barrier_x0, p.barrier_width, p.m, p.hbar);

    case 'trapezoidal'
        T_of_E = @(Ej) analytical_T_trapezoidal(Ej, p.V0, p.barrier_x0, p.barrier_width, p.m, p.hbar);

    case 'double_square'
        T_of_E = @(Ej) analytical_T_double_square(Ej, p.V0, p.barrier_x0, p.barrier_width, p.m, p.hbar);

    case 'gaussian'
        T_of_E = @(Ej) analytical_T_gaussian(Ej, p.V0, p.barrier_x0, p.barrier_width, p.m, p.hbar);

    case 'double_gaussian'
        T_of_E = @(Ej) analytical_T_double_gaussian(Ej, p.V0, p.barrier_x0, p.barrier_width, p.m, p.hbar);

    otherwise
        error('Unknown barrier type: %s', p.barrier_type);

end

E_values = linspace(0.2, 4, 25);
T_sim = zeros(size(E_values));
T_analytical_avg = zeros(size(E_values));

for i = 1:numel(E_values)
    p.E = E_values(i);

    out = simulate_tunneling(p);

    T_sim(i) = out.T;

    T_analytical_avg(i) = analytical_wavepacket_T(p.E, p.sigma, p.m, p.hbar, T_of_E);

    fprintf('E = %5.2f  ->  Numerical T = %.4f  |  Analytical avg T = %.4f\n',p.E, T_sim(i), T_analytical_avg(i));
end

figure('Position', [100 100 900 600]);

plot(E_values, T_sim, 'bo-', 'LineWidth', 1.3, 'MarkerFaceColor', 'b');

hold on;

plot(E_values, T_analytical_avg, 'r--', 'LineWidth', 1.5);

xline(p.V0, 'k:', 'V_0', 'LabelVerticalAlignment', 'bottom');

xlabel('Incident energy E');
ylabel('Transmission coefficient T');

title(sprintf('Quantum tunneling (%s barrier): Numerical vs energy-averaged analytical result', p.barrier_type));

legend('Crank-Nicolson simulation', 'Energy-averaged analytical/transfer-matrix result', 'V_0', 'Location', 'southeast');

grid on;
ylim([0 1.05]);

switch lower(p.barrier_type)

    case 'square'
        exportgraphics(gcf, 'improved_transmission_vs_energy_square.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_square\n');

    case 'triangular'
        exportgraphics(gcf, 'improved_transmission_vs_energy_triangular.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_triangular\n');

    case 'parabolic'
        exportgraphics(gcf, 'improved_transmission_vs_energy_parabolic.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_parabolic\n');

    case 'trapezoidal'
        exportgraphics(gcf, 'improved_transmission_vs_energy_trapezoidal.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_trapezoidal\n');

    case 'double_square'
        exportgraphics(gcf, 'improved_transmission_vs_energy_double_square.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_double_square\n');

    case 'gaussian'
        exportgraphics(gcf, 'improved_transmission_vs_energy_gaussian.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_gaussian\n');

    case 'double_gaussian'
        exportgraphics(gcf, 'improved_transmission_vs_energy_double_gaussian.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png_double_gaussian\n');

    otherwise
        exportgraphics(gcf, 'improved_transmission_vs_energy.png');
        fprintf('\nSaved plot to improved_transmission_vs_energy.png\n');

end


function Tavg = analytical_wavepacket_T(E0, sigma, m, hbar, T_of_E)

    k0 = sqrt(2*m*E0) / hbar;

    sigma_k = 1 / (2*sigma);

    kmin = max(1e-6, k0 - 6*sigma_k);
    kmax = k0 + 6*sigma_k;

    k = linspace(kmin, kmax, 400);

    E = hbar^2 * k.^2 / (2*m);

    Pk = exp(-(k-k0).^2 / (2*sigma_k^2));

    T_E = zeros(size(E));

    for j = 1:numel(E)
        T_E(j) = T_of_E(E(j));
    end

    Tavg = trapz(k, T_E .* Pk) / trapz(k, Pk);

end