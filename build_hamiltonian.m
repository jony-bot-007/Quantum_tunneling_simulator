function H = build_hamiltonian(x, V, hbar, m)
N  = numel(x);
dx = x(2) - x(1);

e  = ones(N, 1);
D2 = spdiags([e, -2*e, e], -1:1, N, N) / dx^2;

Tkin = -(hbar^2 / (2*m)) * D2;    
H = Tkin + spdiags(V(:), 0, N, N);             
end
