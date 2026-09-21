function V = potential_gaussian(x, V0, x1, width)
V = V0 * exp(-4*log(2) * ((x - (x1 + width/2)) / width).^2);
end
