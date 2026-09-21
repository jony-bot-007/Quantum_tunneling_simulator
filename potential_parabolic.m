function V = potential_parabolic(x, V0, x1, width)
xm = x1 + width/2;
V = zeros(size(x));
mask = (x >= x1) & (x <= x1 + width);
V(mask) = V0 * (1 - ((x(mask) - xm)/(width/2)).^2);
end
