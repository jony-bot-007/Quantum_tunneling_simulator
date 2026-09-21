function psi0 = initialize_wavepacket(x, x0, sigma, k0)
norm_const = (2*pi*sigma^2)^(-1/4);
psi0 = norm_const * exp(-(x - x0).^2 / (4*sigma^2)) .* exp(1i * k0 * x);
psi0 = psi0(:); 
end
