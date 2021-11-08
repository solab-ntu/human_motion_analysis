function [] = plot_aligned_data(fname, yy)

if ~contains(fname, "Aligned")
    fname = "Aligned_" + fname;
end

if isfile("Templates\\" + fname + ".mat")
    load("Templates\\" + fname + ".mat");
else
    disp("Can't find " + fname + "!!!");
    return
end

n = length(IK);
l = size(IK(1).data, 1);
t = 0.01*(0:(l-1))';

figure()
subplot(1,2,1)
box on
hold on
for i = 1:n
    plot(t, IK(i).data(:, 1), 'Color', '#0072BD', 'LineWidth', 1.0)
    plot(t, IK(i).data(:, 3), 'Color', '#D95319', 'LineWidth', 1.0)
    plot(t, IK(i).data(:, 5), 'Color', '#EDB120', 'LineWidth', 1.0)
end
hold off
legend("Hip", "Knee", "Ankle")
title(join(split(fname,"_"), 1) + " (Left)")
xlim([t(1), t(end)])
xlabel("Time(sec)")
if nargin == 2
    ylim(yy)
end
ylabel("Angle(degree)")
grid on
grid minor

subplot(1,2,2)
box on
hold on
for i = 1:n
    plot(t, IK(i).data(:, 2), 'Color', '#0072BD', 'LineWidth', 1.0)
    plot(t, IK(i).data(:, 4), 'Color', '#D95319', 'LineWidth', 1.0)
    plot(t, IK(i).data(:, 6), 'Color', '#EDB120', 'LineWidth', 1.0)
end
hold off
legend("Hip", "Knee", "Ankle")
title(join(split(fname,"_"), 1) + " (Right)")
xlim([t(1), t(end)])
xlabel("Time(sec)")
if nargin == 2
    ylim(yy)
end
ylabel("Angle(degree)")
grid on
grid minor

end