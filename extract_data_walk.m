% colunms of data contain:
% ["Time", "Hip R", "Hip L", "Knee R", "Knee L", "Ankle R", "Ankle L"]

dt = 0.01;
save_output = false; % true false
IKR = struct();
IKL = struct();
IKR(1).data = extract_data("Walk_0__01_IK_data", 2.065, 2.657, dt, save_output, '_r');
IKL(1).data = extract_data("Walk_0__01_IK_data", 2.657, 3.265, dt, save_output, '_l');
IKR(2).data = extract_data("Walk_0__02_IK_data", 2.065, 2.632, dt, save_output, '_r');
IKL(2).data = extract_data("Walk_0__02_IK_data", 2.632, 3.240, dt, save_output, '_l');
IKR(3).data = extract_data("Walk_0__03_IK_data", 2.073, 2.632, dt, save_output, '_r');
IKL(3).data = extract_data("Walk_0__03_IK_data", 2.632, 3.182, dt, save_output, '_l');
IKR(4).data = extract_data("Walk_25__01_IK_data", 2.907, 3.557, dt, save_output, '_r');
IKL(4).data = extract_data("Walk_25__01_IK_data", 3.557, 4.207, dt, save_output, '_l');
IKR(5).data = extract_data("Walk_25__02_IK_data", 3.623, 4.315, dt, save_output, '_r');
IKL(5).data = extract_data("Walk_25__02_IK_data", 4.315, 4.957, dt, save_output, '_l');
IKR(6).data = extract_data("Walk_25__03_IK_data", 4.115, 4.748, dt, save_output, '_r');
IKL(6).data = extract_data("Walk_25__03_IK_data", 4.748, 5.332, dt, save_output, '_l');
IKR(7).data = extract_data("Walk_50__01_IK_data", 8.823, 9.398, dt, save_output, '_r');
% IKL(7).data = extract_data("Walk_50__01_IK_data", 9.398, 10.007, dt, save_output);
IKR(8).data = extract_data("Walk_50__02_IK_data", 4.498, 5.098, dt, save_output, '_r');
IKL(7).data = extract_data("Walk_50__02_IK_data", 5.098, 5.707, dt, save_output, '_l');
IKR(9).data = extract_data("Walk_50__03_IK_data", 2.965, 3.615, dt, save_output, '_r');
IKL(8).data = extract_data("Walk_50__03_IK_data", 3.615, 4.273, dt, save_output, '_l');

% add IKL to IKR
% exchange the column (2,3) (4,5) (6,7) for right and left legs
for i = 1:length(IKL)
    IK = IKL(i).data;
    IK(:,[2,3,4,5,6,7]) = IK(:,[3,2,5,4,7,6]);
    IKR(end+1).data = IK;
end

% ["Time", "Hip L", "Hip R", "Knee L", "Knee R", "Ankle L", "Ankle R"]
%%
for i = 1:length(IKR)
    IKR(i).data(:,1) = [];
end

%%

figure()
subplot(2,3,1)
hold on
for i = 1:length(IKR)
    IK = IKR(i).data;
    t = dt*(0:(size(IK,1)-1))';
    plot(t, IK(:, 1), 'Color', '#0072BD', 'LineWidth', 1.0)
end
ylim([-80, 40])
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
    plot(t, IK(:, 3), 'Color', '#D95319', 'LineWidth', 1.0)
end
ylim([-80, 40])
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
    plot(t, IK(:, 5), 'Color', '#EDB120', 'LineWidth', 1.0)
end
ylim([-80, 40])
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
    plot(t, IK(:, 2), 'Color', '#0072BD', 'LineWidth', 1.0)
end
ylim([-80, 40])
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
    plot(t, IK(:, 4), 'Color', '#D95319', 'LineWidth', 1.0)
end
ylim([-80, 40])
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
    plot(t, IK(:, 6), 'Color', '#EDB120', 'LineWidth', 1.0)
end
ylim([-80, 40])
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
box on
grid on

% save("Templates\\Raw_Walk_FL.mat", 'IKR');

%% alignment by offset

% find the median data
min_dist = Inf;
for i = 1:length(IKR)
    d = 0;
    for j = 1:length(IKR)
        if i ~= j
            d = d + myED(IKR(j).data(:, 2:end), IKR(i).data(:, 2:end));
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_offset_avg = IKR(min_id).data(:, 2:end);
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
            [d, IK, wp] = myED(IKR(i).data(:, 2:end), IK_offset_avg);
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
mfPCA(IK_offset)

%% alignment by Dynamic Time Warping

% find the median data

min_dist = Inf;
for i = 1:length(IKR)
    d = 0;
    for j = 1:length(IKR)
        if i ~= j
            d = d + myDTW(IKR(j).data(:, 2:end), IKR(i).data(:, 2:end));
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_DTW_avg = IKR(min_id).data(:, 2:end);

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
        [d, IK, wp] = myDTW(IKR(i).data(:, 2:end), IK_DTW_avg);
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
mfPCA(IK_DTW)

%% alignment by Move-Split-Merge

% find the median data
min_dist = Inf;
for i = 1:length(IKR)
    d = 0;
    for j = 1:length(IKR)
        if i ~= j
            d = d + myMSM(IKR(j).data(:, 2:end), IKR(i).data(:, 2:end));
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_MSM_avg = IKR(min_id).data(:, 2:end);

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
        [d, IK, wp] = myMSM(IKR(i).data(:, 2:end), IK_MSM_avg);
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
mfPCA(IK_MSM)

%% alignment by SRVF + DTW

% f to srvf
IK_srvf = struct();
for i = 1:length(IKR)
    IK_srvf(i).data = diff(IKR(i).data(:,2:end), 1)/dt;
    IK_srvf(i).data = sign(IK_srvf(i).data).*sqrt(abs(IK_srvf(i).data));
end

% find the median data
min_dist = Inf;
for i = 1:length(IK_srvf)
    d = 0;
    for j = 1:length(IK_srvf)
        if i ~= j
            d = d + myDTW(IK_srvf(j).data, IK_srvf(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_srvf_avg = IK_srvf(min_id).data;

% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

pre = Inf;
while true
    IK_srvf_align = struct();
    tmp = 0;
    
    for i = 1:length(IKR)
        [d, IK, wp] = myDTW(IK_srvf(i).data, IK_srvf_avg);
        IK_srvf_align(i).data = IK;
        wp = interp1(linspace(0,1,length(wp)), wp, linspace(0,1,length(wp)+1), 'makima');
        IK_srvf_align(i).wpf = wp;
        tmp = tmp + d;
    end

    IK_srvf_avg = IK_srvf_align(1).data;
    for i = 2:length(IKR)
        IK_srvf_avg = IK_srvf_avg + IK_srvf_align(i).data;
    end
    IK_srvf_avg = IK_srvf_avg/length(IKR);
    
    if abs(tmp - pre)/tmp < 0.01
        break
    else
        pre = tmp;
    end
end

% srvf to f
IK_srvf_rebuild = struct();
for i = 1:length(IK_srvf_align)
    IK_srvf_rebuild(i).data = zeros(size(IK_srvf_align(i).data)+[1, 0]);
    IK_srvf_rebuild(i).data(1,:) = IKR(i).data(1,2:end);
    for j = 1:size(IK_srvf_align(i).data, 1)
        IK_srvf_rebuild(i).data(j+1,:) = IK_srvf_rebuild(i).data(j,:) + IK_srvf_align(i).data(j,:).*abs(IK_srvf_align(i).data(j,:))*dt;
    end
    IK_srvf_rebuild(i).wpf = IK_srvf_align(i).wpf;
end


% IK_srvf_rebuild = struct();
% for i = 1:length(IK_srvf_align)
%     IK_srvf_rebuild(i).data = interp1(linspace(0,1,size(IKR(i).data,1)) , IKR(i).data(:,2:end), IK_srvf_align(i).wpf,'makima');
%     IK_srvf_rebuild(i).wpf = IK_srvf_align(i).wpf;
% end

myPlot(IK_srvf_rebuild, "SRVF + DTW")
var_fs(IK_srvf_rebuild)

%% alignment by SRVF + MSM

% f to srvf
IK_srvf = struct();
for i = 1:length(IKR)
    IK_srvf(i).data = diff(IKR(i).data(:,2:end), 1)/dt;
    IK_srvf(i).data = sign(IK_srvf(i).data).*sqrt(abs(IK_srvf(i).data));
end

% find the median data
min_dist = Inf;
aaa = [];
for i = 1:length(IK_srvf)
    d = 0;
    for j = 1:length(IK_srvf)
        if i ~= j
            d = d + myMSM(IK_srvf(j).data, IK_srvf(i).data);
        end
    end
    if d < min_dist
        min_dist = d;
        min_id = i;
    end
end
IK_srvf_avg = IK_srvf(min_id).data;

% 1. align data by median
% 2. calculate the sum of distane
% 3. compute the mean of aligned data
% 4. align data by mean
% 5. calculate the sum of distane
% 6. if the sum of distane changes less than 1%, then end up. Otherwise, back to 2.

pre = Inf;
for j = 1:2
    IK_srvf_align = struct();
    tmp = 0;
    
    for i = 1:length(IKR)
        [d, IK, wp] = myMSM(IK_srvf(i).data, IK_srvf_avg);
        IK_srvf_align(i).data = IK;
        wp = interp1(linspace(0,1,length(wp)), wp, linspace(0,1,length(wp)+1), 'makima');
        IK_srvf_align(i).wpf = wp;
        tmp = tmp + d;
    end

    IK_srvf_avg = IK_srvf_align(1).data;
    for i = 2:length(IKR)
        IK_srvf_avg = IK_srvf_avg + IK_srvf_align(i).data;
    end
    IK_srvf_avg = IK_srvf_avg/length(IKR);
    
    if abs(tmp - pre)/tmp < 0.01
        break
    else
        pre = tmp;
    end
end

% srvf to f
IK_srvf_rebuild = struct();
for i = 1:length(IK_srvf_align)
    IK_srvf_rebuild(i).data = zeros(size(IK_srvf_align(i).data)+[1, 0]);
    IK_srvf_rebuild(i).data(1,:) = IKR(i).data(1,2:end);
    for j = 1:size(IK_srvf_align(i).data, 1)
        IK_srvf_rebuild(i).data(j+1,:) = IK_srvf_rebuild(i).data(j,:) + IK_srvf_align(i).data(j,:).*abs(IK_srvf_align(i).data(j,:))*dt;
    end
    IK_srvf_rebuild(i).wpf = IK_srvf_align(i).wpf;
end

myPlot(IK_srvf_rebuild, "SRVF + MSM")
var_fs(IK_srvf_rebuild)

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
IK_beta = zeros(size(IKR(i).data, 2)-1, l, length(IKR));
initial_value = zeros(1, size(IKR(i).data, 2)-1); % original mean initial value
for i = 1:length(IKR)
    IK = IKR(i).data(:,2:end);
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
% initial_diff = zeros(1, size(IKR(i).data, 2)-1); % mean initial value of aligned curves
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
mfPCA(IK_ESA)

%%

save_template(IK);
save_aligned_data(IK);
plot_template("Walk_Forward_R", "Walk FR")
plot_template("Walk_Forward_L", "Walk FL")

%% ===================== Aligned Function Plot =========================

function [] = myPlot(IKR, str)

figure()
if str == ""
    subplot(1,2,1)
    box on
    hold on
    for i = 1:length(IKR)
        IK = IKR(i).data;
        plot(IK(:, 1), IK(:, 2), 'Color', '#0072BD', 'LineWidth', 1.0)
        plot(IK(:, 1), IK(:, 4), 'Color', '#D95319', 'LineWidth', 1.0)
        plot(IK(:, 1), IK(:, 6), 'Color', '#EDB120', 'LineWidth', 1.0)
    end
    hold off
    legend("Hip", "Knee", "Ankle")
    title("Original Walk Data (Swing)")
    xlim([0 0.7])
    xlabel("Time(sec)")
    ylim([-80 40])
    ylabel("Angle(degree)")
    grid on
    grid minor

    subplot(1,2,2)
    box on
    hold on
    for i = 1:length(IKR)
        IK = IKR(i).data;
        plot(IK(:, 1), IK(:, 3), 'Color', '#0072BD', 'LineWidth', 1.0)
        plot(IK(:, 1), IK(:, 5), 'Color', '#D95319', 'LineWidth', 1.0)
        plot(IK(:, 1), IK(:, 7), 'Color', '#EDB120', 'LineWidth', 1.0)
    end
    hold off
    legend("Hip", "Knee", "Ankle")
    title("Original Walk Data (Stance)")
    xlim([0 0.7])
    xlabel("Time(sec)")
    ylim([-80 40])
    ylabel("Angle(degree)")
    grid on
    grid minor
    
else
    
    dt = 0.01;
    t = (0:dt:(dt*(size(IKR(1).data,1)-1)))';
    
    subplot(1,2,1)
    box on
    hold on
    for i = 1:length(IKR)
        IK = IKR(i).data;
        plot(t, IK(:, 1), 'Color', '#0072BD', 'LineWidth', 1.0)
        plot(t, IK(:, 3), 'Color', '#D95319', 'LineWidth', 1.0)
        plot(t, IK(:, 5), 'Color', '#EDB120', 'LineWidth', 1.0)
    end
    hold off
    legend("Hip", "Knee", "Ankle")
    title("Walk Data Aligned by " + str + " (Swing)")
    xlim([0, 0.7])
    xlabel("Time(sec)")
    ylim([-80, 40])
    ylabel("Angle(degree)")
    grid on
    grid minor

    subplot(1,2,2)
    box on
    hold on
    for i = 1:length(IKR)
        IK = IKR(i).data;
        plot(t, IK(:, 2), 'Color', '#0072BD', 'LineWidth', 1.0)
        plot(t, IK(:, 4), 'Color', '#D95319', 'LineWidth', 1.0)
        plot(t, IK(:, 6), 'Color', '#EDB120', 'LineWidth', 1.0)
    end
    hold off
    legend("Hip", "Knee", "Ankle")
    title("Walk Data Aligned by " + str + " (Stance)")
    xlim([0, 0.7])
    xlabel("Time(sec)")
    ylim([-80, 40])
    ylabel("Angle(degree)")
    grid on
    grid minor
end
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
    
    save("Templates\\Walk_Forward_R.mat", ...
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
    
    save("Templates\\Walk_Forward_L.mat",... 
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
    
    save("Templates\\Aligned_Walk_Forward_L.mat", 'time', ...
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
    
    save("Templates\\Aligned_Walk_Forward_R.mat", 'time', ...
         'hip_angle_r', 'hip_angle_l', 'knee_angle_r', 'knee_angle_l', 'ankle_angle_r', 'ankle_angle_l');
end

%%

function [] = mfPCA(IK)

    n = length(IK);
    l = size(IK(1).data, 1);
    Data = zeros(l, n, 6);
    mu = zeros(l, 6);
    va = zeros(l, 6);
    
    for i = 1:length(IK)
        Data(:, i, 1) = IK(i).data(:, 1);
        Data(:, i, 2) = IK(i).data(:, 2);
        Data(:, i, 3) = IK(i).data(:, 3);
        Data(:, i, 4) = IK(i).data(:, 4);
        Data(:, i, 5) = IK(i).data(:, 5);
        Data(:, i, 6) = IK(i).data(:, 6);
    end

    mu(:, 1) = mean(Data(:, :, 1), 2);
    mu(:, 2) = mean(Data(:, :, 2), 2);
    mu(:, 3) = mean(Data(:, :, 3), 2);
    mu(:, 4) = mean(Data(:, :, 4), 2);
    mu(:, 5) = mean(Data(:, :, 5), 2);
    mu(:, 6) = mean(Data(:, :, 6), 2);
    
    va(:, 1) = sqrt(mean((Data(:, :, 1) - repmat(mu(:, 1), 1, n)).^2, 2));
    va(:, 2) = sqrt(mean((Data(:, :, 2) - repmat(mu(:, 2), 1, n)).^2, 2));
    va(:, 3) = sqrt(mean((Data(:, :, 3) - repmat(mu(:, 3), 1, n)).^2, 2));
    va(:, 4) = sqrt(mean((Data(:, :, 4) - repmat(mu(:, 4), 1, n)).^2, 2));
    va(:, 5) = sqrt(mean((Data(:, :, 5) - repmat(mu(:, 5), 1, n)).^2, 2));
    va(:, 6) = sqrt(mean((Data(:, :, 6) - repmat(mu(:, 6), 1, n)).^2, 2));
    
    C = [];
    for j = 1:6
       Cj = [];
       for k = 1:6
           Cjk = ((Data(:, :, j) - repmat(mu(:, j), 1, n))./repmat(va(:, j), 1, n))*((Data(:, :, k) - repmat(mu(:, k), 1, n))./repmat(va(:, k), 1, n))';
           Cj = [Cj, Cjk/n];
       end
       C = [C; Cj];
    end
    
    [U, V, ~] = svd(C);
    V = diag(V);
    k = 1;
    while sum(V(1:k)) < sum(V)*0.9
        k = k + 1;
    end
    lambda = V(1:2);
    disp(sum(lambda)/sum(V));
    
    eig_func = zeros(l, 6, k);
    for j = 1:k
        eig_func(:, :, j) = reshape(U(:,j), [l, 6]);
    end   

    sample = zeros(l, 6, 20);

    for i = 1:20
        X = mu;
        for k = 1:length(lambda)
            X = X + normrnd(0, sqrt(lambda(k)))*va.*eig_func(:, :, k);
        end
        sample(:, :, i) = X;
    end

    t = 0.01*(0:(l-1))';

    figure()
    subplot(1, 2, 1)
    box on
    hold on
    for i = 1:20
        plot(t, sample(:, 1, i), 'Color', '#0072BD', 'LineWidth', 1.0)
        plot(t, sample(:, 3, i), 'Color', '#D95319', 'LineWidth', 1.0)
        plot(t, sample(:, 5, i), 'Color', '#EDB120', 'LineWidth', 1.0)
    end
    hold off
    legend("Hip", "Knee", "Ankle")
    title("random sample")
    xlabel("Time(sec)")
    ylabel("Angle(degree)")
    xlim([0, 0.7])
    ylim([-80, 40])
    grid on
    grid minor
    
    subplot(1, 2, 2)
    box on
    hold on
    for i = 1:20
        plot(t, sample(:, 2, i), 'Color', '#0072BD', 'LineWidth', 1.0)
        plot(t, sample(:, 4, i), 'Color', '#D95319', 'LineWidth', 1.0)
        plot(t, sample(:, 6, i), 'Color', '#EDB120', 'LineWidth', 1.0)
    end
    hold off
    legend("Hip", "Knee", "Ankle")
    xlabel("Time(sec)")
    ylabel("Angle(degree)")
    ylim([-80, 40])
    xlim([0, 0.7])
    grid on
    grid minor
    
end