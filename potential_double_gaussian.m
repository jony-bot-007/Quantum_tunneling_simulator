function V = potential_double_gaussian(x, V0, x1, width)
xm = x1 + width/2;
separation = width/4;
sigma = width/10;

x_left  = xm - separation;
x_right = xm + separation;

V = V0 * exp(-((x - x_left).^2)/(2*sigma^2))+ V0 * exp(-((x - x_right).^2)/(2*sigma^2));
end
