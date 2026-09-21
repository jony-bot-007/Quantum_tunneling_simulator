clear all; clc; close all;

p.hbar = 1;
p.m    = 1;      % particle mass (change willingly)

p.xmin = -60;
p.xmax = 60;
p.N    = 800;

p.x0    = -25;            
p.sigma = 2.5;             
p.E     = 1.0;             

p.V0            = 1.5;              
p.barrier_width = 4;                 
p.barrier_x0    = -p.barrier_width/2; 
p.barrier_type  = 'square';     % Any of 'square', 'triangular', 'gaussian', 'parabolic', 'trapezoidal',  'double_square', 'double_gaussian'.

p.dt            = 0.05;
p.nsteps        = 800;
p.store_history = true;

out = simulate_tunneling(p);

fprintf('Results --->\n');
fprintf('Transmission coefficient T = %.4f\n', out.T);
fprintf('Reflection coefficient   R = %.4f\n', out.R);
fprintf('T + R = %.4f  (probability remaining on grid: %.4f)\n', out.T + out.R, out.prob_total);

figure('Position', [100 100 900 700]);
skip   = 5;                                  
Vscale = 0.5 / max(out.V + eps);             

for n = 1:skip:size(out.psi_history, 2)
    psi_n = out.psi_history(:, n);

    subplot(2,1,1);
    plot(out.x, real(psi_n), 'b-', 'LineWidth', 1.2); hold on;
    plot(out.x, imag(psi_n), 'r-', 'LineWidth', 1.2);
    plot(out.x, out.V*Vscale, 'k--', 'LineWidth', 1);
    hold off;
    ylim([-0.6 0.6]);
    xlabel('x'); ylabel('\psi(x,t)');
    legend('Re(\psi)', 'Im(\psi)', 'Barrier (scaled)', 'Location', 'northwest');
    title(sprintf('Wave function  —  t = %.2f', out.t(n)));
    grid on;

    subplot(2,1,2);
    plot(out.x, abs(psi_n).^2, 'm-', 'LineWidth', 1.5); hold on;
    plot(out.x, out.V*Vscale, 'k--', 'LineWidth', 1);
    hold off;
    ylim([0 0.3]);
    xlabel('x'); ylabel('|\psi(x,t)|^2');
    title('Probability density');
    grid on;

    drawnow;
end

exportgraphics(gcf, 'final_state.png');
fprintf('Saved final frame to final_state.png\n');
