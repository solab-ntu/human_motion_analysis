function [] = plot_raw_data(fname, yy)

if ~contains(fname, "Raw")
    fname = "Raw_" + fname;
end

if isfile("Templates\\" + fname + ".mat")
    load("Templates\\" + fname + ".mat");
else
    disp("Can't find " + fname + "!!!");
    return
end

n = length(IK);

figure()
subplot(1,2,1)
box on
hold on
for i = 1:n
    data = IK(i).data;
    t = 0.01*(0:(size(IK(i).data, 1)-1))';
    plot(t, data(:, 1), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, data(:, 3), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, data(:, 5), 'Color', '#EDB120', 'LineWidth', 0.8)
end
hold off
legend("Hip", "Knee", "Ankle")
title(join(split(fname,"_"), 1) + " (Left)")
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
    data = IK(i).data;
    t = 0.01*(0:(size(IK(i).data, 1)-1))';
    plot(t, data(:, 2), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, data(:, 4), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, data(:, 6), 'Color', '#EDB120', 'LineWidth', 0.8)
end
hold off
legend("Hip", "Knee", "Ankle")
title(join(split(fname,"_"), 1) + " (Right)")
xlabel("Time(sec)")
if nargin == 2
    ylim(yy)
end
ylabel("Angle(degree)")
grid on
grid minor

end