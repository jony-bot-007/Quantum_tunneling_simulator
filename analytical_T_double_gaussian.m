function T = analytical_T_double_gaussian(E, V0, x1, width, m, hbar)
T = zeros(size(E));
xm = x1 + width/2;
separation = width/4;
sigma = width/10;
x_left  = xm - separation;
x_right = xm + separation;
xmin = x1 - 8*sigma;
xmax = x1 + width + 8*sigma;
N = 8001;
x = linspace(xmin, xmax, N);
dx = x(2) - x(1);
V = V0 * exp(-(x-x_left).^2/(2*sigma^2)) + V0 * exp(-(x-x_right).^2/(2*sigma^2));

for i = 1:numel(E)
    Ei = E(i);
    if Ei <= 0
        T(i) = 0;
        continue
    end
    k = sqrt(2*m*Ei)/hbar;
    M = eye(2);
    for n = 1:N-1
        Vm = (V(n) + V(n+1))/2;
        if Ei > Vm
            q = sqrt(2*m*(Ei-Vm))/hbar;
            P = [cos(q*dx), sin(q*dx)/q; -q*sin(q*dx), cos(q*dx)];
        elseif Ei < Vm
            kappa = sqrt(2*m*(Vm-Ei))/hbar;
            P = [cosh(kappa*dx), sinh(kappa*dx)/kappa; kappa*sinh(kappa*dx), cosh(kappa*dx)];
        else
            P = [1, dx; 0, 1];
        end
        M = P*M;
    end
    incident = [1; 1i*k];
    reflected = [1; -1i*k];
    transmitted = [1; 1i*k];
    A = [M*reflected, -transmitted];
    b = -M*incident;
    result = A\b;
    r = result(1);
    t = result(2);
    T(i) = abs(t)^2;
end
T(T < 0) = 0;
T(T > 1) = 1;
end