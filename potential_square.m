function V = potential_square(x, V0, x1, width)
V = zeros(size(x));
V(x >= x1 & x <= x1 + width) = V0;
end
