
% colunms of data contain:
% ["lumbar","hip_r","hip_l","knee_r","knee_l","ankle_r","ankle_l"]

dt = 0.01;
time = (0:dt:0.4)';
n = length(time);

IK = struct();
IK(1).data = extract_data("Static_IK_data", 0, Inf, dt);
IK(2).data = extract_data("Deep_squat_0__01_IK_data", 0, 2.4, dt);
IK(3).data = extract_data("Deep_squat_0__01_IK_data", 7.6, Inf, dt);
IK(4).data = extract_data("Deep_squat_0__02_IK_data", 0, 3.6, dt);
IK(5).data = extract_data("Deep_squat_0__02_IK_data", 8.8, Inf, dt);
IK(6).data = extract_data("Deep_squat_0__03_IK_data", 8.3, Inf, dt);
IK(7).data = extract_data("Deep_squat_25__01_IK_data", 7.5, 8.8, dt);
IK(8).data = extract_data("Deep_squat_25__01_IK_data", 13.7, Inf, dt);
IK(9).data = extract_data("Deep_squat_25__02_IK_data", 2.6, 3.1, dt);
IK(10).data = extract_data("Deep_squat_25__02_IK_data", 7.2, 7.8, dt);
IK(11).data = extract_data("Deep_squat_25__03_IK_data", 2.7, 3.1, dt);
IK(12).data = extract_data("Deep_squat_25__03_IK_data", 7.0, 7.6, dt);
IK(13).data = extract_data("Deep_squat_50__01_IK_data", 7.6, 8.2, dt);
IK(14).data = extract_data("Deep_squat_50__02_IK_data", 8.3, 8.8, dt);
IK(15).data = extract_data("Deep_squat_50__02_IK_data", 13.1, Inf, dt);
IK(16).data = extract_data("Deep_squat_50__03_IK_data", 7.6, 8.3, dt);
IK(17).data = extract_data("Stair_L_h_a1_IK_data", 0, 1.8, dt);
IK(18).data = extract_data("Stair_L_h_a2_IK_data", 0, 1.8, dt);
IK(19).data = extract_data("Stair_L_h_a3_IK_data", 0, 1.5, dt);
IK(20).data = extract_data("Stair_L_h_a4_IK_data", 0, 2, dt);
IK(21).data = extract_data("Stair_L_h_a5_IK_data", 0, 1.8, dt);
IK(22).data = extract_data("Stair_L_n_a1_IK_data", 0, 2, dt);
IK(23).data = extract_data("Stair_L_n_a2_IK_data", 0, 1.4, dt);
IK(24).data = extract_data("Stair_L_n_a3_IK_data", 0, 2, dt);
IK(25).data = extract_data("Stair_L_n_a4_IK_data", 0, 1.5, dt);
IK(26).data = extract_data("Stair_L_n_a5_IK_data", 0, 1.6, dt);
IK(27).data = extract_data("Stair_R_h_a1_IK_data", 0, 1.4, dt);
IK(28).data = extract_data("Stair_R_h_a2_IK_data", 0, 2, dt);
IK(29).data = extract_data("Stair_R_h_a3_IK_data", 0, 1.5, dt);
IK(30).data = extract_data("Stair_R_h_a4_IK_data", 0.6, 1.8, dt);
IK(31).data = extract_data("Stair_R_h_a5_IK_data", 0, 2, dt);
IK(32).data = extract_data("Stair_R_n_a1_IK_data", 0, 1.6, dt);
IK(33).data = extract_data("Stair_R_n_a2_IK_data", 0, 1.4, dt);
IK(34).data = extract_data("Stair_R_n_a3_IK_data", 0, 2, dt);
IK(35).data = extract_data("Stair_R_n_a4_IK_data", 0, 1.4, dt);
IK(36).data = extract_data("Stair_R_n_a5_IK_data", 0, 1.4, dt);
IK(37).data = extract_data("Stair_L_h_d1_IK_data", 0, 1.8, dt);
IK(38).data = extract_data("Stair_L_h_d2_IK_data", 0, 1.7, dt);
IK(39).data = extract_data("Stair_L_h_d3_IK_data", 0, 1.6, dt);
IK(40).data = extract_data("Stair_L_h_d4_IK_data", 0, 1.8, dt);
IK(41).data = extract_data("Stair_L_h_d5_IK_data", 0, 1.7, dt);
IK(42).data = extract_data("Stair_L_n_d1_IK_data", 0, 2.1, dt);
IK(43).data = extract_data("Stair_L_n_d2_IK_data", 0, 1.3, dt);
IK(44).data = extract_data("Stair_L_n_d3_IK_data", 0, 1.6, dt);
IK(45).data = extract_data("Stair_L_n_d4_IK_data", 0, 2.1, dt);
IK(46).data = extract_data("Stair_L_n_d5_IK_data", 0, 2.2, dt);
IK(47).data = extract_data("Stair_R_h_d1_IK_data", 0, 1.8, dt);
IK(48).data = extract_data("Stair_R_h_d2_IK_data", 0, 1.9, dt);
IK(49).data = extract_data("Stair_R_h_d3_IK_data", 0, 1.9, dt);
IK(50).data = extract_data("Stair_R_h_d4_IK_data", 0, 1.7, dt);
IK(51).data = extract_data("Stair_R_h_d5_IK_data", 0, 1.7, dt);
IK(52).data = extract_data("Stair_R_n_d1_IK_data", 0, 1.9, dt);
IK(53).data = extract_data("Stair_R_n_d2_IK_data", 0, 1.6, dt);
IK(54).data = extract_data("Stair_R_n_d3_IK_data", 0, 1.9, dt);
IK(55).data = extract_data("Stair_R_n_d4_IK_data", 0, 1.6, dt);
IK(56).data = extract_data("Stair_R_n_d5_IK_data", 0, 1.7, dt);

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

%%

while length(ankle_angle) < 17835
    hip_angle = [hip_angle; repmat(hip_angle(end), 17835 - length(hip_angle), 1)];
    knee_angle = [knee_angle; repmat(knee_angle(end), 17835 - length(knee_angle), 1)];
    ankle_angle = [ankle_angle; repmat(ankle_angle(end), 17835 - length(ankle_angle), 1)];
end

hip_angle = reshape(hip_angle, 41, []);
knee_angle = reshape(knee_angle, 41, []);
ankle_angle = reshape(ankle_angle, 41, []);

% hip_angle = repmat(mean(hip_angle, 1), n, 1);
% knee_angle = repmat(mean(knee_angle, 1), n, 1);
% ankle_angle = repmat(mean(ankle_angle, 1), n, 1);

%%

hip_angle_mean = mean(hip_angle, 2);
knee_angle_mean = mean(knee_angle, 2);
ankle_angle_mean= mean(ankle_angle, 2);

hip_angle_std = std(hip_angle, 0, 2);
knee_angle_std = std(knee_angle, 0, 2);
ankle_angle_std = std(ankle_angle, 0, 2);


%%

r = randi([1, size(hip_angle, 2)],1,50);

figure()
subplot(1,3,1)
hold on
plot(time, hip_angle(:,r), 'Color', '#0072BD', 'LineWidth', 0.8);
ylim([-30, 30])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(1,3,2)
hold on
plot(time, knee_angle(:,r), 'Color', '#D95319', 'LineWidth', 0.8);
ylim([-30, 30])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(1,3,3)
hold on
plot(time, ankle_angle(:,r), 'Color', '#EDB120', 'LineWidth', 0.8);
ylim([-30, 30])
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

save("Templates\\Stand.mat", ...
     'time', 'hip_angle_r', 'hip_angle_r_std', 'hip_angle_l', 'hip_angle_l_std', ...
     'knee_angle_r', 'knee_angle_r_std', 'knee_angle_l', 'knee_angle_l_std', ...
     'ankle_angle_r', 'ankle_angle_r_std', 'ankle_angle_l', 'ankle_angle_l_std');

plot_template("Stand")

%%

hip_angle_r = hip_angle;
hip_angle_l = hip_angle;
knee_angle_r = knee_angle;
knee_angle_l = knee_angle;
ankle_angle_r = ankle_angle;
ankle_angle_l = ankle_angle;

save("Templates\\Aligned_Stand.mat", 'time', ...
     'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');

plot_aligned_data("Aligned_Stand")

%%

hip_angle_r = hip_angle;
hip_angle_l = hip_angle;
knee_angle_r = knee_angle;
knee_angle_l = knee_angle;
ankle_angle_r = ankle_angle;
ankle_angle_l = ankle_angle;

save("Templates\\Raw_Stand.mat", 'time', ...
     'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');
