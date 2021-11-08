% colunms of data contain:
% ["Time", "Hip R", "Hip L", "Knee R", "Knee L", "Ankle R", "Ankle L"]

dt = 0.01;
save_output = false; % true false
IKD = struct();
IKU = struct();
IKD(1).data = extract_data("Deep_squat_0__01_IK_data", 2.5, 4.8, dt);
IKD(2).data = extract_data("Deep_squat_0__02_IK_data", 3.6, 5.8, dt);
IKD(3).data = extract_data("Deep_squat_0__03_IK_data", 3.7, 5.6, dt);
IKD(4).data = extract_data("Deep_squat_25__01_IK_data", 3.1, 5.5, dt);
IKD(5).data = extract_data("Deep_squat_25__01_IK_data", 8.9, 10.8, dt);
IKD(6).data = extract_data("Deep_squat_25__02_IK_data", 3.2, 4.8, dt);
IKD(7).data = extract_data("Deep_squat_25__02_IK_data", 7.8, 10, dt);
IKD(8).data = extract_data("Deep_squat_25__03_IK_data", 3.2, 4.8, dt);
IKD(9).data = extract_data("Deep_squat_25__03_IK_data", 7.6, 9.8, dt);
IKD(10).data = extract_data("Deep_squat_50__01_IK_data", 8.2, 10.4, dt);
IKD(11).data = extract_data("Deep_squat_50__02_IK_data", 4.6, 6.2, dt);
IKD(12).data = extract_data("Deep_squat_50__02_IK_data", 8.8, 10.8, dt);
IKD(13).data = extract_data("Deep_squat_50__03_IK_data", 8.2, 10.4, dt);

IKU(1).data = extract_data("Deep_squat_0__01_IK_data", 5.6, 7.6, dt);
IKU(2).data = extract_data("Deep_squat_0__02_IK_data", 6.9, 8.8, dt);
IKU(3).data = extract_data("Deep_squat_0__03_IK_data", 6.5, 8.4, dt);
IKU(4).data = extract_data("Deep_squat_25__01_IK_data", 5.5, 7.5, dt);
IKU(5).data = extract_data("Deep_squat_25__01_IK_data", 11.8, 13.6, dt);
IKU(6).data = extract_data("Deep_squat_25__02_IK_data", 5.4, 7.1, dt);
IKU(7).data = extract_data("Deep_squat_25__02_IK_data", 10.4, 12, dt);
IKU(8).data = extract_data("Deep_squat_25__03_IK_data", 5.3, 7, dt);
IKU(9).data = extract_data("Deep_squat_25__03_IK_data", 10.2, 12, dt);
IKU(10).data = extract_data("Deep_squat_50__01_IK_data", 10.8, 12.3, dt);
IKU(11).data = extract_data("Deep_squat_50__02_IK_data", 6.7, 8.3, dt);
IKU(12).data = extract_data("Deep_squat_50__02_IK_data", 11.7, 13.1, dt);
IKU(13).data = extract_data("Deep_squat_50__03_IK_data", 5.8, 7.5, dt);
IKU(14).data = extract_data("Deep_squat_50__03_IK_data", 10.8, 12.6, dt);

% add IKL to IKD
% exchange the column (2,3) (4,5) (6,7) for right and left legs
for i = 1:length(IKU)
    IK = IKU(i).data;
    IK(:,[2,3,4,5,6,7]) = flip(IK(:,[3,2,5,4,7,6]),1);
    IKD(end+1).data = IK;
end

l = length(IKD);
for i = 1:l
    IKD(i).data(:,1) = [];
    IKD(end+1).data = IKD(i).data(:,[2,4,6]);
    IKD(i).data(:,[2,4,6]) = [];
end

%%

% myPlot(IKD, "")

figure()
subplot(1,3,1)
hold on
for i = 1:length(IKD)
    IK = IKD(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 1), 'Color', '#0072BD', 'LineWidth', 1.0)
end
ylim([-140, 120])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(1,3,2)
hold on
for i = 1:length(IKD)
    IK = IKD(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 2), 'Color', '#D95319', 'LineWidth', 1.0)
end
ylim([-140, 120])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

subplot(1,3,3)
hold on
for i = 1:length(IKD)
    IK = IKD(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 3), 'Color', '#EDB120', 'LineWidth', 1.0)
end
ylim([-140, 120])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

% save("Templates\\Raw_Squat_Down.mat", 'IKD');

%% alignment by offset

% find the median data
min_dist = Inf;
for i = 1:length(IKD)
    d = 0;
    for j = 1:length(IKD)
        if i ~= j
            d = d + myED(IKD(j).data, IKD(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_offset_avg = IKD(min_id).data;
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
    
    for i = 1:length(IKD)
        if i == id
            IK_offset(i).data = IK_offset_avg;
            id = 0;
        else
            [d, IK, wp] = myED(IKD(i).data, IK_offset_avg);
            IK_offset(i).data = IK;
            IK_offset(i).wpf = wp;
            tmp = tmp + d;
        end
    end

    IK_offset_avg = IK_offset(1).data;
    for i = 2:length(IKD)
        IK_offset_avg = IK_offset_avg + IK_offset(i).data;
    end
    IK_offset_avg = IK_offset_avg/length(IKD);
    
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
for i = 1:length(IKD)
    d = 0;
    for j = 1:length(IKD)
        if i ~= j
            d = d + myDTW(IKD(j).data, IKD(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_DTW_avg = IKD(min_id).data;

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
    
    for i = 1:length(IKD)
        [d, IK, wp] = myDTW(IKD(i).data, IK_DTW_avg);
        IK_DTW(i).data = IK;
        IK_DTW(i).wpf = wp;
        tmp = tmp + d;
    end

    IK_DTW_avg = IK_DTW(1).data;
    for i = 2:length(IKD)
        IK_DTW_avg = IK_DTW_avg + IK_DTW(i).data;
    end
    IK_DTW_avg = IK_DTW_avg/length(IKD);
    
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
for i = 1:length(IKD)
    % i
    d = 0;
    for j = 1:length(IKD)
        if i ~= j
            d = d + myMSM(IKD(j).data, IKD(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_MSM_avg = IKD(min_id).data;

% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

% pre = Inf;
for j = 1:1
    IK_MSM = struct();
    tmp = 0;
    
    for i = 1:length(IKD)
        [d, IK, wp] = myMSM(IKD(i).data, IK_MSM_avg);
        IK_MSM(i).data = IK;
        IK_MSM(i).wpf = wp;
        tmp = tmp + d;
        % i 
    end

    IK_MSM_avg = IK_MSM(1).data;
    for i = 2:length(IKD)
        IK_MSM_avg = IK_MSM_avg + IK_MSM(i).data;
    end
    IK_MSM_avg = IK_MSM_avg/length(IKD);
    
%     if abs(tmp - pre)/tmp < 0.01
%         break
%     else
%         pre = tmp;
%     end
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
for i = 1:length(IKD)
    l = l + size(IKD(i).data, 1);
end
l = round(l/length(IKD));

% initialize input curve
IK_beta = zeros(size(IKD(i).data, 2), l, length(IKD));
initial_value = zeros(1, size(IKD(i).data, 2)); % original mean initial value
for i = 1:length(IKD)
    IK = IKD(i).data;
    IK = interp1(linspace(0,1,size(IKD(i).data, 1))', IK, linspace(0,1,l)', 'makima');
    IK_beta(:,:,i) = IK';
    initial_value = initial_value + IK(1,:);
end
initial_value = initial_value/length(IKD);

% computing aligned curves
ESA_obj = fdacurve(IK_beta, false, l);
ESA_objn = ESA_obj.karcher_mean();

% % output aligned curves
% IK_ESA = struct();
% initial_diff = zeros(1, size(IKD(i).data, 2)); % mean initial value of aligned curves
% for i = 1:length(IKD)
%     IK_ESA(i).data = (ESA_objn.betan(:,:,i))';
%     IK_ESA(i).wpf = ESA_objn.gams(:,i);
%     initial_diff = initial_diff + IK_ESA(i).data(1,:);
% end
% initial_diff = initial_value - initial_diff/length(IKD); % differece of initial value
% 
% IK_ESA_avg = (ESA_objn.beta_mean)' + ones(l,1)*initial_diff;
% for i = 1:length(IKD)
%     IK_ESA(i).data = IK_ESA(i).data + ones(l,1)*initial_diff;
% end

IK_ESA = struct();
initial_diff = zeros(1, size(IKR(i).data, 2)); % mean initial value of aligned curves
for i = 1:length(IKD)
    IK_ESA(i).wpf = ESA_objn.gams(:,i);
    IK_ESA(i).data = interp1(linspace(0,1,size(IK_beta(:,:,i),2)),IK_beta(:,:,i)',IK_ESA(i).wpf,'makima');
end

myPlot(IK_ESA, "ESA")
var_fs(IK_ESA)

%%

save_pattern(IK)
save_aligned_data(IK)

%%

save_pattern(IK_ESA)
save_aligned_data(IK_ESA)
plot_template("Squat_Down", "Squat Down Template")
plot_template("Squat_Up", "Squat Up Template")

%% ===================== Aligned Function Plot =========================

function [] = myPlot(IKD, str)

dt = 0.01;

figure()
box on
hold on
for i = 1:length(IKD)
    IK = IKD(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 1), 'Color', '#0072BD', 'LineWidth', 1.0)
    plot(t, IK(:, 2), 'Color', '#D95319', 'LineWidth', 1.0)
    plot(t, IK(:, 3), 'Color', '#EDB120', 'LineWidth', 1.0)
end
hold off
legend("Hip", "Knee", "Ankle")

if str == ""
    title("Original Squat Data")
else
    title("Squat Data Aligned by " + str)
end

xlim([0, 2.5])
xlabel("Time(sec)")
ylim([-140, 120])
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

function [] = save_pattern(IK)
    
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
    
    dt = 0.01;
    time = dt*(0:(size(IKm,1)-1))';
    
    % save mat
%     hip_angle_r = IKm(:,1);
%     hip_angle_r_std = IKs(:,1);
%     hip_angle_l = IKm(:,1);
%     hip_angle_l_std = IKs(:,1);
%     knee_angle_r = IKm(:,2);
%     knee_angle_r_std = IKs(:,2);
%     knee_angle_l = IKm(:,2);
%     knee_angle_l_std = IKs(:,2);
%     ankle_angle_r = IKm(:,3);
%     ankle_angle_r_std = IKs(:,3);
%     ankle_angle_l = IKm(:,3);
%     ankle_angle_l_std = IKs(:,3);
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
    
    save("Templates\\Squat_Down.mat", ...
         'time', 'hip_angle_r', 'hip_angle_r_std', 'hip_angle_l', 'hip_angle_l_std', ...
         'knee_angle_r', 'knee_angle_r_std', 'knee_angle_l', 'knee_angle_l_std', ...
         'ankle_angle_r', 'ankle_angle_r_std', 'ankle_angle_l', 'ankle_angle_l_std');
    
    hip_angle_r = flip(hip_angle_r);
    hip_angle_r_std = flip(hip_angle_r_std);
    hip_angle_l = flip(hip_angle_l);
    hip_angle_l_std = flip(hip_angle_l_std);
    knee_angle_r = flip(knee_angle_r);
    knee_angle_r_std = flip(knee_angle_r_std);
    knee_angle_l = flip(knee_angle_l);
    knee_angle_l_std = flip(knee_angle_l_std);
    ankle_angle_r = flip(ankle_angle_r);
    ankle_angle_r_std = flip(ankle_angle_r_std);
    ankle_angle_l = flip(ankle_angle_l);
    ankle_angle_l_std = flip(ankle_angle_l_std);
    
    save("Templates\\Squat_Up.mat",... 
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
%         hip_angle_r = [hip_angle_r, IK(i).data(:,1)];
%         hip_angle_l = [hip_angle_l, IK(i).data(:,1)];
%         knee_angle_r = [knee_angle_r, IK(i).data(:,2)];
%         knee_angle_l = [knee_angle_l, IK(i).data(:,2)];
%         ankle_angle_r = [ankle_angle_r, IK(i).data(:,3)];
%         ankle_angle_l = [ankle_angle_l, IK(i).data(:,3)];
        hip_angle_r = [hip_angle_r, IK(i).data(:,1)];
        hip_angle_l = [hip_angle_l, IK(i).data(:,2)];
        knee_angle_r = [knee_angle_r, IK(i).data(:,3)];
        knee_angle_l = [knee_angle_l, IK(i).data(:,4)];
        ankle_angle_r = [ankle_angle_r, IK(i).data(:,5)];
        ankle_angle_l = [ankle_angle_l, IK(i).data(:,6)];
    end
    
    save("Templates\\Aligned_Squat_Down.mat", 'time', ...
         'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');

    hip_angle_r = flip(hip_angle_r, 1);
    hip_angle_l = flip(hip_angle_l, 1);
    knee_angle_r = flip(knee_angle_r, 1);
    knee_angle_l = flip(knee_angle_l, 1);
    ankle_angle_r = flip(ankle_angle_r, 1);
    ankle_angle_l = flip(ankle_angle_l, 1);
    
    save("Templates\\Aligned_Squat_Up.mat", 'time', ...
         'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');
end
