function out = simulate_tunneling(p)
x  = linspace(p.xmin, p.xmax, p.N)';
dx = x(2) - x(1);

V = build_potential(x, p);

k0  = sqrt(2*p.m*p.E) / p.hbar;
psi = initialize_wavepacket(x, p.x0, p.sigma, k0);
psi = psi / sqrt(sum(abs(psi).^2) * dx);

H = build_hamiltonian(x, V, p.hbar, p.m);
N = p.N;
I = speye(N);

A = I + 1i * p.dt / (2*p.hbar) * H;
B = I - 1i * p.dt / (2*p.hbar) * H;

dA = decomposition(A);   % factorize once, reuse every step

if p.store_history
    psi_history = zeros(N, p.nsteps + 1);
    psi_history(:, 1) = psi;
end
t = (0:p.nsteps) * p.dt;

for n = 1:p.nsteps
    psi = crank_nicolson_step(psi, dA, B);
    psi([1 end]) = 0;   % psi=0 at edges

    if p.store_history
        psi_history(:, n+1) = psi;
    end
end

barrier_x1 = p.barrier_x0 + p.barrier_width;
[T, R, prob_total] = compute_TR(x, psi, p.barrier_x0, barrier_x1);

out.x          = x;
out.t          = t;
out.V          = V;
out.psi_final  = psi;
out.T          = T;
out.R          = R;
out.prob_total = prob_total;

if p.store_history
    out.psi_history = psi_history;
end
end
