
x = 0:0.1:2*pi;
y = sin(x);

yy = exp(x/3) - 1;
yy = yy/max(yy);

yr = makima(x, y, yy*2*pi);

figure()
plot(x, y, 'LineWidth', 1.2);
xlim([0 2*pi])

figure()
plot(x, yr, 'LineWidth', 1.2);
xlim([0 2*pi])

figure()
plot(x/max(x), yy, 'Color', '#D95319', 'LineWidth', 1.2);
xlim([0 1])
ylim([0 1])
