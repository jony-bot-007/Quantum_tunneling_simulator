function V = build_potential(x, p)
switch lower(p.barrier_type)
    case 'square'
        V = potential_square(x, p.V0, p.barrier_x0, p.barrier_width);
    case 'triangular'
        V = potential_triangular(x, p.V0, p.barrier_x0, p.barrier_width);
    case 'gaussian'
        V = potential_gaussian(x, p.V0, p.barrier_x0, p.barrier_width);
    case 'parabolic'
        V = potential_parabolic(x, p.V0, p.barrier_x0, p.barrier_width);
    case 'trapezoidal'
        V = potential_trapezoidal(x, p.V0, p.barrier_x0, p.barrier_width);
    case 'double_square'
        V = potential_double_square(x, p.V0, p.barrier_x0, p.barrier_width);
    case 'double_gaussian'
        V = potential_double_gaussian(x, p.V0, p.barrier_x0, p.barrier_width);
    otherwise
        error('build_potential:badBarrier', ...
            ['Unknown barrier type "%s". ''Use square, triangular, gaussian, parabolic, ''trapezoidal, double_square, or double_gaussian.'], p.barrier_type);
end
end
