function V = potential_double_square(x, V0, x1, width)
V = zeros(size(x));

barrier_width = 0.25 * width;
gap_width     = 0.50 * width;

b1_start = x1;
b1_end   = x1 + barrier_width;
b2_start = x1 + barrier_width + gap_width;
b2_end   = x1 + width;

mask1 = (x >= b1_start) & (x <= b1_end);
mask2 = (x >= b2_start) & (x <= b2_end);

V(mask1) = V0;
V(mask2) = V0;
end
