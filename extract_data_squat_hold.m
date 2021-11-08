% colunms of data contain:
% ["lumbar","hip_r","hip_l","knee_r","knee_l","ankle_r","ankle_l"]
% 
% PLEASE CHECK Line 40 in extract_data.m

dt = 0.01;
time = (0:dt:0.4)';
n = length(time);

IK = struct();
IK(1).data = extract_data("Deep_squat_0__01_IK_data", 4.8, 5.6, dt);
IK(2).data = extract_data("Deep_squat_0__02_IK_data", 5.9, 6.9, dt);
IK(3).data = extract_data("Deep_squat_0__03_IK_data", 5.6, 6.6, dt);
IK(4).data = extract_data("Deep_squat_25__01_IK_data", 5.4, 5.6, dt);
IK(5).data = extract_data("Deep_squat_25__01_IK_data", 10.8, 12, dt);
IK(6).data = extract_data("Deep_squat_25__02_IK_data", 4.8, 5.6, dt);
IK(7).data = extract_data("Deep_squat_25__02_IK_data", 9.9, 10.6, dt);
IK(8).data = extract_data("Deep_squat_25__03_IK_data", 4.8, 5.4, dt);
IK(9).data = extract_data("Deep_squat_25__03_IK_data", 9.6, 10.4, dt);
IK(10).data = extract_data("Deep_squat_50__01_IK_data", 10, 11, dt);
IK(11).data = extract_data("Deep_squat_50__02_IK_data", 6, 7, dt);
IK(12).data = extract_data("Deep_squat_50__02_IK_data", 10.6, 11.9, dt);
IK(13).data = extract_data("Deep_squat_50__03_IK_data", 5.3, 6, dt);
IK(14).data = extract_data("Deep_squat_50__03_IK_data", 10, 11, dt);

%%

hip_angle_r = [];
hip_angle_l = [];
knee_angle_r = [];
knee_angle_l = [];
ankle_angle_r = [];
ankle_angle_l = [];

for i = 1:length(IK)
    hip_angle_r = [hip_angle_r, repmat(mean(IK(i).data(:,2)), n, 1)];
    hip_angle_l = [hip_angle_l, repmat(mean(IK(i).data(:,3)), n, 1)];
    knee_angle_r = [knee_angle_r, repmat(mean(IK(i).data(:,4)), n, 1)];
    knee_angle_l = [knee_angle_l, repmat(mean(IK(i).data(:,5)), n, 1)];
    ankle_angle_r = [ankle_angle_r, repmat(mean(IK(i).data(:,6)), n, 1)];
    ankle_angle_l = [ankle_angle_l, repmat(mean(IK(i).data(:,7)), n, 1)];
end

hip_angle = [hip_angle_r, hip_angle_l];
knee_angle = [knee_angle_r, knee_angle_l];
ankle_angle = [ankle_angle_r, ankle_angle_l];


hip_angle_mean = mean(hip_angle, 2);
knee_angle_mean = mean(knee_angle, 2);
ankle_angle_mean = mean(ankle_angle, 2);

hip_angle_std = std(hip_angle, 0, 2);
knee_angle_std = std(knee_angle, 0, 2);
ankle_angle_std = std(ankle_angle, 0, 2);

%%

hip_angle_r = [];
hip_angle_l = [];
knee_angle_r = [];
knee_angle_l = [];
ankle_angle_r = [];
ankle_angle_l = [];

for i = 1:length(IK)
    hip_angle_r = [hip_angle_r; IK(i).data(:,2)];
    hip_angle_l = [hip_angle_l; IK(i).data(:,3)];
    knee_angle_r = [knee_angle_r; IK(i).data(:,4)];
    knee_angle_l = [knee_angle_l; IK(i).data(:,5)];
    ankle_angle_r = [ankle_angle_r; IK(i).data(:,6)];
    ankle_angle_l = [ankle_angle_l; IK(i).data(:,7)];
end

hip_angle = [hip_angle_r; hip_angle_l];
knee_angle = [knee_angle_r; knee_angle_l];
ankle_angle = [ankle_angle_r; ankle_angle_l];

while length(ankle_angle) < 2460
    hip_angle = [hip_angle; repmat(hip_angle(end), 2460 - length(hip_angle), 1)];
    knee_angle = [knee_angle; repmat(knee_angle(end), 2460 - length(knee_angle), 1)];
    ankle_angle = [ankle_angle; repmat(ankle_angle(end), 2460 - length(ankle_angle), 1)];
end

hip_angle = reshape(hip_angle, 41, []);
knee_angle = reshape(knee_angle, 41, []);
ankle_angle = reshape(ankle_angle, 41, []);

%%

figure()
subplot(1,3,1)
hold on
plot(time, hip_angle, 'Color', '#0072BD', 'LineWidth', 0.8);
ylim([-140, 120])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(1,3,2)
hold on
plot(time, knee_angle, 'Color', '#D95319', 'LineWidth', 0.8);
ylim([-140, 120])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(1,3,3)
hold on
plot(time, ankle_angle, 'Color', '#EDB120', 'LineWidth', 0.8);
ylim([-140, 120])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

%%

hip_angle_r = hip_angle_mean;
hip_angle_l = hip_angle_mean;
knee_angle_r = knee_angle_mean;
knee_angle_l = knee_angle_mean;
ankle_angle_r = ankle_angle_mean;
ankle_angle_l = ankle_angle_mean;

hip_angle_r_std = hip_angle_std;
hip_angle_l_std = hip_angle_std;
knee_angle_r_std = knee_angle_std;
knee_angle_l_std = knee_angle_std;
ankle_angle_r_std = ankle_angle_std;
ankle_angle_l_std = ankle_angle_std;
     
save("Templates\\Squat_hold.mat", ...
     'time', 'hip_angle_r', 'hip_angle_r_std', 'hip_angle_l', 'hip_angle_l_std', ...
     'knee_angle_r', 'knee_angle_r_std', 'knee_angle_l', 'knee_angle_l_std', ...
     'ankle_angle_r', 'ankle_angle_r_std', 'ankle_angle_l', 'ankle_angle_l_std');

plot_template("Squat_hold")
 
%%

hip_angle_r = hip_angle;
hip_angle_l = hip_angle;
knee_angle_r = knee_angle;
knee_angle_l = knee_angle;
ankle_angle_r = ankle_angle;
ankle_angle_l = ankle_angle;

save("Templates\\Aligned_Squat_Hold.mat", 'time', ...
     'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');

 plot_aligned_data("Aligned_Squat_Hold")

 %%
 
hip_angle_r = hip_angle;
hip_angle_l = hip_angle;
knee_angle_r = knee_angle;
knee_angle_l = knee_angle;
ankle_angle_r = ankle_angle;
ankle_angle_l = ankle_angle;

save("Templates\\Raw_Squat_Hold.mat", 'time', ...
     'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');

 plot_raw_data("Raw_Squat_Hold")
 