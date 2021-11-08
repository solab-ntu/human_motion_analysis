% colunms of data contain:
% ["lumbar","hip_r","hip_l","knee_r","knee_l","ankle_r","ankle_l"]
% 
% PLEASE CHECK Line 40 in extract_data.m

dt = 0.01;
IKR = struct();
IKL = struct();

IKR(1).data = extract_data("Stair_L_n_d1_IK_data", 2.88, 3.62, dt);
IKR(2).data = extract_data("Stair_L_n_d1_IK_data", 4.25, 4.85, dt);
IKR(3).data = extract_data("Stair_L_n_d2_IK_data", 2.35, 3.04, dt);
IKR(4).data = extract_data("Stair_L_n_d2_IK_data", 3.66, 4.21, dt);
IKR(5).data = extract_data("Stair_L_n_d3_IK_data", 2.76, 3.40, dt);
IKR(6).data = extract_data("Stair_L_n_d3_IK_data", 4.02, 4.57, dt);
IKR(7).data = extract_data("Stair_L_n_d4_IK_data", 3.48, 4.13, dt);
IKR(8).data = extract_data("Stair_L_n_d4_IK_data", 4.74, 5.31, dt);
IKR(9).data = extract_data("Stair_L_n_d5_IK_data", 3.26, 3.95, dt);
IKR(10).data = extract_data("Stair_L_n_d5_IK_data", 4.53, 5.07, dt);
IKR(11).data = extract_data("Stair_L_h_d1_IK_data", 3.02, 3.74, dt);
IKR(12).data = extract_data("Stair_L_h_d1_IK_data", 4.40, 4.98, dt);
IKR(13).data = extract_data("Stair_L_h_d2_IK_data", 3.18, 3.90, dt);
IKR(14).data = extract_data("Stair_L_h_d2_IK_data", 4.50, 5.07, dt);
IKR(15).data = extract_data("Stair_L_h_d3_IK_data", 2.91, 3.60, dt);
IKR(16).data = extract_data("Stair_L_h_d3_IK_data", 4.21, 4.77, dt);
IKR(17).data = extract_data("Stair_L_h_d4_IK_data", 3.13, 3.80, dt);
IKR(18).data = extract_data("Stair_L_h_d4_IK_data", 4.41, 5.00, dt);
IKR(19).data = extract_data("Stair_L_h_d5_IK_data", 2.90, 3.66, dt);
IKR(20).data = extract_data("Stair_L_h_d5_IK_data", 4.30, 4.86, dt);
IKR(21).data = extract_data("Stair_R_n_d1_IK_data", 3.69, 4.28, dt);
IKR(22).data = extract_data("Stair_R_n_d2_IK_data", 3.41, 4.02, dt);
IKR(23).data = extract_data("Stair_R_n_d3_IK_data", 3.61, 4.21, dt);
IKR(24).data = extract_data("Stair_R_n_d4_IK_data", 3.27, 3.84, dt);
IKR(25).data = extract_data("Stair_R_n_d5_IK_data", 3.40, 4.02, dt);
IKR(26).data = extract_data("Stair_R_h_d1_IK_data", 3.76, 4.42, dt);
IKR(27).data = extract_data("Stair_R_h_d2_IK_data", 3.76, 4.46, dt);
IKR(28).data = extract_data("Stair_R_h_d3_IK_data", 3.84, 4.53, dt);
IKR(29).data = extract_data("Stair_R_h_d4_IK_data", 3.62, 4.25, dt);
IKR(30).data = extract_data("Stair_R_h_d5_IK_data", 3.61, 4.21, dt);

IKL(1).data = extract_data("Stair_L_n_d1_IK_data", 3.62, 4.25, dt);
IKL(2).data = extract_data("Stair_L_n_d2_IK_data", 3.04, 3.66, dt);
IKL(3).data = extract_data("Stair_L_n_d3_IK_data", 3.40, 4.02, dt);
IKL(4).data = extract_data("Stair_L_n_d4_IK_data", 4.13, 4.74, dt);
IKL(5).data = extract_data("Stair_L_n_d5_IK_data", 3.95, 4.53, dt);
IKL(6).data = extract_data("Stair_L_h_d1_IK_data", 3.74, 4.40, dt);
IKL(7).data = extract_data("Stair_L_h_d2_IK_data", 3.90, 4.50, dt);
IKL(8).data = extract_data("Stair_L_h_d3_IK_data", 3.60, 4.21, dt);
IKL(9).data = extract_data("Stair_L_h_d4_IK_data", 3.80, 4.41, dt);
IKL(10).data = extract_data("Stair_L_h_d5_IK_data", 3.66, 4.30, dt);
IKL(11).data = extract_data("Stair_R_n_d1_IK_data", 3.02, 3.69, dt);
IKL(12).data = extract_data("Stair_R_n_d1_IK_data", 4.28, 4.81, dt);
IKL(13).data = extract_data("Stair_R_n_d2_IK_data", 2.75, 3.41, dt);
IKL(14).data = extract_data("Stair_R_n_d2_IK_data", 4.02, 4.56, dt);
IKL(15).data = extract_data("Stair_R_n_d3_IK_data", 2.96, 3.61, dt);
IKL(16).data = extract_data("Stair_R_n_d3_IK_data", 4.21, 4.75, dt);
IKL(17).data = extract_data("Stair_R_n_d4_IK_data", 2.67, 3.27, dt);
IKL(18).data = extract_data("Stair_R_n_d4_IK_data", 3.84, 4.37, dt);
IKL(19).data = extract_data("Stair_R_n_d5_IK_data", 2.73, 3.40, dt);
IKL(20).data = extract_data("Stair_R_n_d5_IK_data", 4.02, 4.61, dt);
IKL(21).data = extract_data("Stair_R_h_d1_IK_data", 3.02, 3.76, dt);
IKL(22).data = extract_data("Stair_R_h_d1_IK_data", 4.42, 5.03, dt);
IKL(23).data = extract_data("Stair_R_h_d2_IK_data", 3.05, 3.76, dt);
IKL(24).data = extract_data("Stair_R_h_d2_IK_data", 4.46, 5.05, dt);
IKL(25).data = extract_data("Stair_R_h_d3_IK_data", 3.11, 3.84, dt);
IKL(26).data = extract_data("Stair_R_h_d3_IK_data", 4.53, 5.07, dt);
IKL(27).data = extract_data("Stair_R_h_d4_IK_data", 2.95, 3.62, dt);
IKL(28).data = extract_data("Stair_R_h_d4_IK_data", 4.25, 4.81, dt);
IKL(29).data = extract_data("Stair_R_h_d5_IK_data", 2.96, 3.61, dt);
IKL(30).data = extract_data("Stair_R_h_d5_IK_data", 4.21, 4.71, dt);

% add IKL to IKR
% exchange the column (2,3) (4,5) (6,7) for right and left legs
for i = 1:length(IKL)
    IK = IKL(i).data;
    IK(:,[2,3,4,5,6,7]) = IK(:,[3,2,5,4,7,6]);
    IKR(end+1).data = IK();
end

for i = 1:length(IKR)
    IKR(i).data(:,1) = [];
end

%%

% myPlot(IKR, "")

figure()
subplot(2,3,1)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 1), 'Color', '#0072BD', 'LineWidth', 0.8)
end
ylim([-120, 80])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(2,3,2)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 3), 'Color', '#D95319', 'LineWidth', 0.8)
end
ylim([-120, 80])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(2,3,3)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 5), 'Color', '#EDB120', 'LineWidth', 0.8)
end
ylim([-120, 80])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(2,3,4)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 2), 'Color', '#0072BD', 'LineWidth', 0.8)
end
ylim([-120, 80])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(2,3,5)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 4), 'Color', '#D95319', 'LineWidth', 0.8)
end
ylim([-120, 80])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(2,3,6)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 6), 'Color', '#EDB120', 'LineWidth', 0.8)
end
ylim([-120, 80])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

save("Templates\\Raw_Stair_DL.mat", 'IKR');

%% alignment by offset

% find the median data
min_dist = Inf;
for i = 1:length(IKR)
    d = 0;
    for j = 1:length(IKR)
        if i ~= j
            d = d + myED(IKR(j).data, IKR(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_offset_avg = IKR(min_id).data;
id = min_id;

% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

pre = Inf;
while true
    clear IK_offset
    IK_offset = struct();
    tmp = 0;
    
    for i = 1:length(IKR)
        if i == id
            IK_offset(i).data = IK_offset_avg;
            id = 0;
        else
            [d, IK, wp] = myED(IKR(i).data, IK_offset_avg);
            IK_offset(i).data = IK;
            IK_offset(i).wpf = wp;
            tmp = tmp + d;
        end
    end

    IK_offset_avg = IK_offset(1).data;
    for i = 2:length(IKR)
        IK_offset_avg = IK_offset_avg + IK_offset(i).data;
    end
    IK_offset_avg = IK_offset_avg/length(IKR);
    
    if abs(tmp - pre)/tmp < 0.01
        break
    else
        pre = tmp;
    end
end


myPlot(IK_offset, "Offset")
var_fs(IK_offset)

%% alignment by Dynamic Time Warping

% find the median data

min_dist = Inf;
for i = 1:length(IKR)
    d = 0;
    for j = 1:length(IKR)
        if i ~= j
            d = d + myDTW(IKR(j).data, IKR(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_DTW_avg = IKR(min_id).data;

% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

pre = Inf;
while true
    IK_DTW = struct();
    tmp = 0;
    
    for i = 1:length(IKR)
        [d, IK, wp] = myDTW(IKR(i).data, IK_DTW_avg);
        IK_DTW(i).data = IK;
        IK_DTW(i).wpf = wp;
        tmp = tmp + d;
    end

    IK_DTW_avg = IK_DTW(1).data;
    for i = 2:length(IKR)
        IK_DTW_avg = IK_DTW_avg + IK_DTW(i).data;
    end
    IK_DTW_avg = IK_DTW_avg/length(IKR);
    
    if abs(tmp - pre)/tmp < 0.01
        break
    else
        pre = tmp;
    end
end

myPlot(IK_DTW, "DTW")
var_fs(IK_DTW)

%% alignment by Move-Split-Merge

% find the median data
min_dist = Inf;
for i = 1:length(IKR)
    d = 0;
    for j = 1:length(IKR)
        if i ~= j
            d = d + myMSM(IKR(j).data, IKR(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_MSM_avg = IKR(min_id).data;

% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

pre = Inf;
for j = 1:1
    IK_MSM = struct();
    tmp = 0;
    
    for i = 1:length(IKR)
        [d, IK, wp] = myMSM(IKR(i).data, IK_MSM_avg);
        IK_MSM(i).data = IK;
        IK_MSM(i).wpf = wp;
        tmp = tmp + d;
    end

    IK_MSM_avg = IK_MSM(1).data;
    for i = 2:length(IKR)
        IK_MSM_avg = IK_MSM_avg + IK_MSM(i).data;
    end
    IK_MSM_avg = IK_MSM_avg/length(IKR);
    
    if abs(tmp - pre)/tmp < 0.01
        break
    else
        pre = tmp;
    end
end

myPlot(IK_MSM, "MSM")
var_fs(IK_MSM)

%% alignment by Elastic Shape Analysis
% using tool fdasrvf developed by J. Derek Tucker, refering to arXiv:1212.1791v2
% fdacurve class provides alignment methods for curves in R^n using the SRVF framework
%
% obj = fdacurve(beta, closed, N, scale)
% input:
%   beta       (n,T,K) matrix defining n dimensional and T samples of K curves
%   closed     true or false if closed curve
%   N          resample curve to N points
%   scale      include scale (optional: false by defualt)
% output:
%   obj        a fdacurve object
%
% fdacurve class
% properties:
%   beta       (n,T,K) matrix of curve
%   q          (n,T,K) matrix of srvfs
%   betan      aligned curves
%   qn         aligned srvfs
%   beta_mean  karcher mean curve
%   q_mean     karcher mean srvf
%   gams       warping functions
% methods:
%   obj = karcher_mean(obj, option)

% find IK_beta for same number of samples

% mean length of curves
l = 0;
for i = 1:length(IKR)
    l = l + size(IKR(i).data, 1);
end
l = round(l/length(IKR));

% initialize input curve
IK_beta = zeros(size(IKR(i).data, 2), l, length(IKR));
initial_value = zeros(1, size(IKR(i).data, 2)); % original mean initial value
for i = 1:length(IKR)
    IK = IKR(i).data;
    IK = interp1(linspace(0,1,size(IKR(i).data, 1))', IK, linspace(0,1,l)', 'makima');
    IK_beta(:,:,i) = IK';
    initial_value = initial_value + IK(1,:);
end
initial_value = initial_value/length(IKR);

% computing aligned curves
ESA_obj = fdacurve(IK_beta, false, l);
ESA_objn = ESA_obj.karcher_mean();

% % output aligned curves
% IK_ESA = struct();
% initial_diff = zeros(1, size(IKR(i).data, 2)); % mean initial value of aligned curves
% for i = 1:length(IKR)
%     IK_ESA(i).data = (ESA_objn.betan(:,:,i))';
%     IK_ESA(i).wpf = ESA_objn.gams(:,i);
%     initial_diff = initial_diff + IK_ESA(i).data(1,:);
% end
% initial_diff = initial_value - initial_diff/length(IKR); % differece of initial value
% 
% IK_ESA_avg = (ESA_objn.beta_mean)' + ones(l,1)*initial_diff;
% for i = 1:length(IKR)
%     IK_ESA(i).data = IK_ESA(i).data + ones(l,1)*initial_diff;
% end

IK_ESA = struct();
initial_diff = zeros(1, size(IKR(i).data, 2)); % mean initial value of aligned curves
for i = 1:length(IKR)
    IK_ESA(i).wpf = ESA_objn.gams(:,i);
    IK_ESA(i).data = interp1(linspace(0,1,size(IK_beta(:,:,i),2)),IK_beta(:,:,i)',IK_ESA(i).wpf,'makima');
end

myPlot(IK_ESA, "ESA")
var_fs(IK_ESA)

%%

save_template(IK);
save_aligned_data(IK);
plot_template("Stair_Descent_R", "Stair DR Template")
plot_template("Stair_Descent_L", "Stair DL Template")


%% ===================== Aligned Function Plot =========================

function [] = myPlot(IKR, str)

dt = 0.01;

figure()
subplot(1,2,1)
box on
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 1), 'Color', '#0072BD', 'LineWidth', 1.0)
    plot(t, IK(:, 3), 'Color', '#D95319', 'LineWidth', 1.0)
    plot(t, IK(:, 5), 'Color', '#EDB120', 'LineWidth', 1.0)
end
hold off
legend("Hip", "Knee", "Ankle")
if str == ""
    title("Original Stair Descent Data (Swing)")
else
    title("Stair Descent Data Aligned by " + str + " (Swing)")
end
xlim([0, 0.8])
xlabel("Time(sec)")
ylim([-120, 60])
ylabel("Angle(degree)")
grid on
grid minor

subplot(1,2,2)
box on
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 2), 'Color', '#0072BD', 'LineWidth', 1.0)
    plot(t, IK(:, 4), 'Color', '#D95319', 'LineWidth', 1.0)
    plot(t, IK(:, 6), 'Color', '#EDB120', 'LineWidth', 1.0)
end
hold off
legend("Hip", "Knee", "Ankle")
if str == ""
    title("Original Stair Descent Data (Stance)")
else
    title("Stair Descent Data Aligned by " + str + " (Stance)")
end
xlim([0, 0.8])
xlabel("Time(sec)")
ylim([-120, 60])
ylabel("Angle(degree)")
grid on
grid minor

end

function v = var_fs(IK)
% vy   amplitube variation
% xy   amplitube variation

    IKmw = struct();
    IKxm = zeros(size(IK(1).data));
    IKym = IKxm;
    l = 0;
    for i = 1:length(IK)
        IKym = IKym + IK(i).data;
        l = l + size(IKym,1);
    end
    IKym = IKym/length(IK);
    l = l/length(IK);
    
    vy = 0;
    for i = 1:length(IK)
        IKmw(i).data = interp1(linspace(0,1,size(IKym,1))', IKym, IK(i).wpf, 'makima');
        IKxm = IKxm + IKmw(i).data;
        vy = vy + sum(sum((IK(i).data - IKym).^2));
    end
    vy = vy/length(IK);
    
    IKxm = IKxm/length(IK);
    
    vx = 0;
    for i = 1:length(IK)
        vx = vx + sum(sum((IKmw(i).data - IKxm).^2));
    end
    vx = vx/length(IK);
    
    v = [vy, vx]/l;
end

function [] = save_template(IK)
    
    IKm = IK(1).data;
    for i =2:length(IK)
        IKm = IKm + IK(i).data;
    end
    IKm = IKm/length(IK);
    
    IKs = zeros(size(IKm));
    for i =2:length(IK)
        IKs = IKs + (IK(i).data - IKm).^2;
    end
    IKs = sqrt(IKs/(length(IK)-1));
    
    a = 1;
    dt = 0.01;
    t = 0:dt:(dt*(size(IKm,1)-1)); 
    
    % save mat
    time = t';
    hip_angle_r = IKm(:,1);
    hip_angle_r_std = IKs(:,1);
    hip_angle_l = IKm(:,2);
    hip_angle_l_std = IKs(:,2);
    knee_angle_r = IKm(:,3);
    knee_angle_r_std = IKs(:,3);
    knee_angle_l = IKm(:,4);
    knee_angle_l_std = IKs(:,4);
    ankle_angle_r = IKm(:,5);
    ankle_angle_r_std = IKs(:,5);
    ankle_angle_l = IKm(:,6);
    ankle_angle_l_std = IKs(:,6);
    
    save("Templates\\Stair_Descent_R.mat", ...
         'time', 'hip_angle_r', 'hip_angle_r_std', 'hip_angle_l', 'hip_angle_l_std', ...
         'knee_angle_r', 'knee_angle_r_std', 'knee_angle_l', 'knee_angle_l_std', ...
         'ankle_angle_r', 'ankle_angle_r_std', 'ankle_angle_l', 'ankle_angle_l_std');
    
    hip_angle_r = IKm(:,2);
    hip_angle_r_std = IKs(:,2);
    hip_angle_l = IKm(:,1);
    hip_angle_l_std = IKs(:,1);
    knee_angle_r = IKm(:,4);
    knee_angle_r_std = IKs(:,4);
    knee_angle_l = IKm(:,3);
    knee_angle_l_std = IKs(:,3);
    ankle_angle_r = IKm(:,6);
    ankle_angle_r_std = IKs(:,6);
    ankle_angle_l = IKm(:,5);
    ankle_angle_l_std = IKs(:,5);
    
    save("Templates\\Stair_Descent_L.mat",... 
         'time', 'hip_angle_r', 'hip_angle_r_std', 'hip_angle_l', 'hip_angle_l_std', ...
         'knee_angle_r', 'knee_angle_r_std', 'knee_angle_l', 'knee_angle_l_std', ...
         'ankle_angle_r', 'ankle_angle_r_std', 'ankle_angle_l', 'ankle_angle_l_std');
end

function [] = save_aligned_data(IK)
    
    dt = 0.01;
    time = dt*(0:(size(IK(1).data,1)-1))';

    hip_angle_r = [];
    hip_angle_l = [];
    knee_angle_r = [];
    knee_angle_l = [];
    ankle_angle_r = [];
    ankle_angle_l = [];

    for i =2:length(IK)
        hip_angle_r = [hip_angle_r, IK(i).data(:,1)];
        hip_angle_l = [hip_angle_l, IK(i).data(:,2)];
        knee_angle_r = [knee_angle_r, IK(i).data(:,3)];
        knee_angle_l = [knee_angle_l, IK(i).data(:,4)];
        ankle_angle_r = [ankle_angle_r, IK(i).data(:,5)];
        ankle_angle_l = [ankle_angle_l, IK(i).data(:,6)];
    end
    
    save("Templates\\Aligned_Stair_Descent_R.mat", 'time', ...
         'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');

    hip_angle_r = hip_angle_r + hip_angle_l;
    knee_angle_r = knee_angle_r + knee_angle_l;
    ankle_angle_r = ankle_angle_r + ankle_angle_l;
    
    hip_angle_l = hip_angle_r - hip_angle_l;
    knee_angle_l = knee_angle_r - knee_angle_l;
    ankle_angle_l = ankle_angle_r - ankle_angle_l;
    
    hip_angle_r = hip_angle_r - hip_angle_l;
    knee_angle_r = knee_angle_r - knee_angle_l;
    ankle_angle_r = ankle_angle_r - ankle_angle_l;
    
    save("Templates\\Aligned_Stair_Descent_L.mat", 'time', ...
         'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');
end