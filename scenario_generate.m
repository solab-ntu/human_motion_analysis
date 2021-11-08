
template_names = ["Stand";
                  "Squat_Hold";
                  "Squat_Down";
                  "Squat_Up";
                  "Walk_FL";
                  "Walk_FR";
                  "Stair_AL";
                  "Stair_AR";
                  "Stair_DL";
                  "Stair_DR"];

id_series = [1, 1, 1, 3, 2, 4, 1, 1, 5, 6, 5, 6, 1, 1, 8, 7, 8, 7, 1, 1, 10, 9, 10, 9, 1];
temp = load_mfpca(template_names);

scenario = [];
for i = 1:length(id_series)
    motion = motion_gen(temp(id_series(i)));
    
    if i > 1
        t = [1:5, 16:20]';
        tq = (1:20)';
        motion_s = [scenario((end - 4):end, :); motion(1:5, :)];
        motion_s = interp1(t, motion_s, tq, 'spline');
        scenario = [scenario; motion_s(6:15, :); motion];
    else
        scenario = [scenario; motion];
    end
end

%%

t = 0.01*((1:size(scenario, 1))'-1);

figure()
% clf

box on
hold on
plot(t, scenario(:, 1), 'Color', '#0072BD', 'LineWidth', 1.2)
plot(t, scenario(:, 2), 'Color', '#0072BD', 'LineStyle', '--', 'LineWidth', 1.0)
plot(t, scenario(:, 3), 'Color', '#D95319', 'LineWidth', 1.2)
plot(t, scenario(:, 4), 'Color', '#D95319', 'LineStyle', '--', 'LineWidth', 1.0)
plot(t, scenario(:, 5), 'Color', '#EDB120', 'LineWidth', 1.2)
plot(t, scenario(:, 6), 'Color', '#EDB120', 'LineStyle', '--', 'LineWidth', 1.0)
hold off
legend("Hip L", "Hip R", "Knee L", "Knee R", "Ankle L", "Ankle R")
xlim([t(1) t(end)])
ylim([-150 120])
xlabel("Time(sec)")
ylabel("Angle(degree)")
xticks(t(1):2:t(end))
grid on

% save("Scenario\\Scenario_2.mat", 'scenario');

%%

scenario = load_angle_motion("Scenario\\Scenario_old.mat");
scenario = interp1(scenario(:, 1), scenario(:, 2:7), scenario(1, 1):0.01:scenario(end, 1), 'makima');

%%

load("Scenario\\Scenario_1.mat");

scenario(1800:2101, :) = [];
scenario(1150:1250, :) = [];
scenario(1:100, :) = [];

motion = motion_gen(temp(id_series(1)));
t = [1:5, 16:20]';
tq = (1:20)';
motion_s = [scenario((end - 4):end, :); motion(1:5, :)];
motion_s = interp1(t, motion_s, tq, 'spline');
scenario = [scenario; motion_s(6:15, :); motion];

%%

load("Scenario\\Scenario_1.mat")
l = size(scenario, 1);
t = 0.01*((1:l)'-1);

scenario = scenario + normrnd(0, 1, [l, 6]);

% y2 = scenario + normrnd(0, 1, [l, 6]);
% y3 = scenario + normrnd(0, 2, [l, 6]);
% y4 = scenario + normrnd(0, 3, [l, 6]);

% figure()
% subplot(4,1,1)
% plot(t, scenario)
% subplot(4,1,2)
% plot(t, y2)
% subplot(4,1,3)
% plot(t, y3)
% subplot(4,1,4)
% plot(t, y4)

% figure()
clf

box on
hold on
plot(t, scenario(:, 1), 'Color', '#0072BD', 'LineWidth', 1.2)
plot(t, scenario(:, 2), 'Color', '#0072BD', 'LineStyle', '--', 'LineWidth', 1.0)
plot(t, scenario(:, 3), 'Color', '#D95319', 'LineWidth', 1.2)
plot(t, scenario(:, 4), 'Color', '#D95319', 'LineStyle', '--', 'LineWidth', 1.0)
plot(t, scenario(:, 5), 'Color', '#EDB120', 'LineWidth', 1.2)
plot(t, scenario(:, 6), 'Color', '#EDB120', 'LineStyle', '--', 'LineWidth', 1.0)
hold off
legend("Hip L", "Hip R", "Knee L", "Knee R", "Ankle L", "Ankle R")
xlim([t(1) t(end)])
ylim([-150 120])
xlabel("Time(sec)")
ylabel("Angle(degree)")
xticks(t(1):2:t(end))
grid on

%%

x = 0:0.02:(2*pi);
l = length(x);

y1 = sin(x);
y2 = y1 + normrnd(0, 1, [1, l]);
y3 = y1 + movmean(normrnd(0, 1, [1, l]), 5);

figure()
subplot(3,1,1)
plot(x, y1)
subplot(3,1,2)
plot(x, y2)
subplot(3,1,3)
plot(x, y3)

%% motion generate example

l = size(temp(3).mu, 1);
n = 30;
motion = struct();
gamma = zeros(l, n);
for i = 1:n
    [motion(i).data, g]= motion_gen(temp(3));
    gamma(:, i) = g;
end

figure()
subplot(1,2,1)
hold on
for i = 1:n
    t = 0:1/(l-1):1;
    plot(t, gamma(:, i), 'LineWidth', 0.8)
end
axis([0 1 0 1])
hold off
grid on
box on

subplot(1,2,2)
hold on
for i = 1:n
    Data = motion(i).data;
    t = 0.01*((1:size(Data, 1))' -1);
    plot(t, Data(:, 1), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data(:, 3), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data(:, 5), 'Color', '#EDB120', 'LineWidth', 0.8)
end
legend("Hip", "Knee", "Ankle")
xlabel("time(sec)")
ylabel("angle(degree)")
hold off
grid on
box on

%%

function [motion, gamma] = motion_gen(temp)

    mu = temp.mu;
    std = temp.std;
    eig_func = temp.eig_func;
    lambda = temp.lambda;
    l = size(mu, 1);
    
    score = normrnd(0, sqrt(lambda)/2);
    motion = zeros(l, 6);
    
    for i = 1:length(score)
        motion = motion + score(i)*eig_func(:, :, i);
    end
    
    motion = motion.*std + mu;
    
    a = normrnd(0, 0.05);
    b = normrnd(0, 0.03);
    c = normrnd(1.07, 0.08);
    c = round(l*c);
    
    t = (0:1/(l-1):1)';
    tq = (0:1/(c-1):1)';
    gamma = tq + a*sin(tq*pi) + a*sin(1*tq*pi);
    motion = interp1(t, motion, gamma, 'makima');
end

function temp = load_mfpca(template_names)
    
    temp = struct();
    
    for i = 1:size(template_names, 1)
        
        fin = load("Templates\\Aligned_" + template_names(i) + ".mat");
        n = length(fin.IK); % number of samples
        l = size(fin.IK(1).data, 1);
        
        Data = zeros(l, n, 6);
        mu = zeros(l, 6);
        std = zeros(l, 6);
        for j = 1:n
            Data(:, j, 1) = fin.IK(j).data(:, 1);
            Data(:, j, 2) = fin.IK(j).data(:, 2);
            Data(:, j, 3) = fin.IK(j).data(:, 3);
            Data(:, j, 4) = fin.IK(j).data(:, 4);
            Data(:, j, 5) = fin.IK(j).data(:, 5);
            Data(:, j, 6) = fin.IK(j).data(:, 6);
        end
        mu(:, 1) = mean(Data(:, j, 1), 2);
        mu(:, 2) = mean(Data(:, j, 2), 2);
        mu(:, 3) = mean(Data(:, j, 3), 2);
        mu(:, 4) = mean(Data(:, j, 4), 2);
        mu(:, 5) = mean(Data(:, j, 5), 2);
        mu(:, 6) = mean(Data(:, j, 6), 2);
        std(:, 1) = sqrt(mean((Data(:, :, 1) - repmat(mu(:, 1), 1, n)).^2, 2));
        std(:, 2) = sqrt(mean((Data(:, :, 2) - repmat(mu(:, 2), 1, n)).^2, 2));
        std(:, 3) = sqrt(mean((Data(:, :, 3) - repmat(mu(:, 3), 1, n)).^2, 2));
        std(:, 4) = sqrt(mean((Data(:, :, 4) - repmat(mu(:, 4), 1, n)).^2, 2));
        std(:, 5) = sqrt(mean((Data(:, :, 5) - repmat(mu(:, 5), 1, n)).^2, 2));
        std(:, 6) = sqrt(mean((Data(:, :, 6) - repmat(mu(:, 6), 1, n)).^2, 2));
        
        C = [];
        for j = 1:6
           Cj = [];
           for k = 1:6
               Cjk = ((Data(:, :, j) - repmat(mu(:, j), 1, n))./repmat(std(:, j), 1, n))*((Data(:, :, k) - repmat(mu(:, k), 1, n))./repmat(std(:, k), 1, n))';
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

        eig_func = zeros(l, 6, k);
        for j = 1:k
            eig_func(:, :, j) = reshape(U(:,j), [l, 6]);
        end
        
        temp(i).mu = mu;
        temp(i).std = std;
        temp(i).eig_func = eig_func;
        temp(i).lambda = V(1:k);
    end
end
