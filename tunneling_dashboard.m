function tunneling_dashboard

fig = uifigure('Name', 'Quantum Tunneling Simulator', 'Position', [80 80 1100 680]);

mainGL = uigridlayout(fig, [1 2]);
mainGL.ColumnWidth = {260, '1x'};

ctrlGL = uigridlayout(mainGL, [10 2]);
ctrlGL.Layout.Row = 1; ctrlGL.Layout.Column = 1;
ctrlGL.RowHeight = repmat({35}, 1, 9);
ctrlGL.RowHeight{10} = '1x';

uilabel(ctrlGL, 'Text', 'Barrier type', 'FontWeight', 'bold');
ddType = uidropdown(ctrlGL, ...
    'Items', { ...
    'square', ...
    'triangular', ...
    'gaussian', ...
    'parabolic', ...
    'trapezoidal', ...
    'double_square', ...
    'double_gaussian'});

uilabel(ctrlGL, 'Text', 'Barrier height V0');
efV0 = uieditfield(ctrlGL, 'numeric', 'Value', 1.5, 'Limits', [0 10]);

uilabel(ctrlGL, 'Text', 'Barrier width L');
efL = uieditfield(ctrlGL, 'numeric', 'Value', 4, 'Limits', [0.5 20]);

uilabel(ctrlGL, 'Text', 'Particle energy E');
efE = uieditfield(ctrlGL, 'numeric', 'Value', 1.0, 'Limits', [0.05 10]);

uilabel(ctrlGL, 'Text', 'Particle mass m');
efM = uieditfield(ctrlGL, 'numeric', 'Value', 1.0, 'Limits', [0.1 10]);

uilabel(ctrlGL, 'Text', 'Wave packet width \sigma');
efSigma = uieditfield(ctrlGL, 'numeric', 'Value', 2.5, 'Limits', [0.5 8]);

btnRun = uibutton(ctrlGL, 'Text', 'Run Simulation', 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @runSim);
btnRun.Layout.Column = [1 2];

lblResult = uilabel(ctrlGL, 'Text', 'T = --    R = --');
lblResult.Layout.Column = [1 2];
lblResult.FontWeight = 'bold';

lblCheck = uilabel(ctrlGL, 'Text', '');
lblCheck.Layout.Column = [1 2];

filler = uilabel(ctrlGL, 'Text', '');
filler.Layout.Row = 10; filler.Layout.Column = [1 2];

plotGL = uigridlayout(mainGL, [2 1]);
plotGL.Layout.Row = 1; plotGL.Layout.Column = 2;

axWave = uiaxes(plotGL);
axWave.Layout.Row = 1;
title(axWave, 'Wave function \psi(x,t)');
xlabel(axWave, 'x'); ylabel(axWave, '\psi');

axProb = uiaxes(plotGL);
axProb.Layout.Row = 2;
title(axProb, 'Probability density |\psi(x,t)|^2');
xlabel(axProb, 'x'); ylabel(axProb, '|\psi|^2');

    function runSim(~, ~)
        p.hbar = 1;
        p.m    = efM.Value;

        p.xmin = -60; p.xmax = 60; p.N = 700;

        p.x0    = -25;
        p.sigma = efSigma.Value;
        p.E     = efE.Value;

        p.V0            = efV0.Value;
        p.barrier_width = efL.Value;
        p.barrier_x0    = -efL.Value / 2;
        p.barrier_type  = ddType.Value;

        p.dt            = 0.05;
        p.nsteps        = 700;
        p.store_history = true;

        btnRun.Enable = 'off';
        lblResult.Text = 'Running...';
        lblCheck.Text  = '';
        drawnow;

        out = simulate_tunneling(p);

        Vscale = 0.5 / max(out.V + eps);
        skip   = 8;

        for n = 1:skip:size(out.psi_history, 2)
            psi_n = out.psi_history(:, n);

            cla(axWave);
            plot(axWave, out.x, real(psi_n), 'b-', 'LineWidth', 1.2); hold(axWave, 'on');
            plot(axWave, out.x, imag(psi_n), 'r-', 'LineWidth', 1.2);
            plot(axWave, out.x, out.V*Vscale, 'k--');
            hold(axWave, 'off');
            ylim(axWave, [-0.6 0.6]);
            legend(axWave, {'Re(\psi)', 'Im(\psi)', 'Barrier'}, 'Location', 'northwest');

            cla(axProb);
            plot(axProb, out.x, abs(psi_n).^2, 'm-', 'LineWidth', 1.5); hold(axProb, 'on');
            plot(axProb, out.x, out.V*Vscale, 'k--');
            hold(axProb, 'off');
            ylim(axProb, [0 0.3]);

            drawnow;
        end

        lblResult.Text = sprintf('T = %.4f    R = %.4f', out.T, out.R);
        lblCheck.Text  = sprintf('T + R = %.4f  (prob. conserved: %.4f)', ...
            out.T + out.R, out.prob_total);

        btnRun.Enable = 'on';
    end
end
