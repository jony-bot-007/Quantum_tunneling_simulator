function V = potential_trapezoidal(x, V0, x1, width)
V = zeros(size(x));
xm1 = x1 + 0.25*width;
xm2 = x1 + 0.75*width;

leftMask = (x >= x1) & (x < xm1);
middleMask = (x >= xm1) & (x <= xm2);
rightMask = (x > xm2) & (x <= x1 + width);

V(leftMask) = V0 *(x(leftMask) - x1) / (0.25*width);
V(middleMask) = V0;
V(rightMask) = V0 *(x1 + width - x(rightMask)) / (0.25*width);
end
