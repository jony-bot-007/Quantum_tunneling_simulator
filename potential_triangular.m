function V = potential_triangular(x, V0, x1, width)
V = zeros(size(x));
xm = x1 + width/2;

leftMask  = (x >= x1) & (x < xm);
rightMask = (x >= xm) & (x <= x1 + width);

V(leftMask)  = V0 * (x(leftMask) - x1) / (width/2);
V(rightMask) = V0 * (x1 + width - x(rightMask)) / (width/2);
end
