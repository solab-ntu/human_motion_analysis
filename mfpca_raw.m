
% load("Templates\\Raw_Stand.mat");
% load("Templates\\Raw_Squat_Hold.mat");
% load("Templates\\Raw_Squat_Down.mat");
% load("Templates\\Raw_Walk_FL.mat");
% load("Templates\\Raw_Stair_AL.mat");
% load("Templates\\Raw_Stair_DL.mat");

% load("Templates\\Aligned_Stand.mat");
% load("Templates\\Aligned_Squat_Hold.mat");
% load("Templates\\Aligned_Squat_Down.mat");
% load("Templates\\Aligned_Walk_FL.mat");
% load("Templates\\Aligned_Stair_AL.mat");
load("Templates\\Aligned_Stair_DL.mat");

%%

n = length(IK); % number of samples
l = 0;
for i = 1:n
    l = l + size(IK(i).data, 1);
end
l = round(l/n); % mean length of samples

Data = zeros(l, n, 6);
mu = zeros(l, 6);
std = zeros(l, 6);
for i = 1:n
    DataIn = IK(i).data;
    t = 0:1/(size(DataIn, 1) - 1):1;
    tq = 0:1/(l - 1):1;
    DataIn = interp1(t', DataIn, tq', 'makima');
    Data(:, i, 1) = DataIn(:, 1);
    Data(:, i, 2) = DataIn(:, 2);
    Data(:, i, 3) = DataIn(:, 3);
    Data(:, i, 4) = DataIn(:, 4);
    Data(:, i, 5) = DataIn(:, 5);
    Data(:, i, 6) = DataIn(:, 6);
end

mu(:, 1) = mean(Data(:, :, 1), 2);
mu(:, 2) = mean(Data(:, :, 2), 2);
mu(:, 3) = mean(Data(:, :, 3), 2);
mu(:, 4) = mean(Data(:, :, 4), 2);
mu(:, 5) = mean(Data(:, :, 5), 2);
mu(:, 6) = mean(Data(:, :, 6), 2);
std(:, 1) = sqrt(mean((Data(:, :, 1) - repmat(mu(:, 1), 1, n)).^2, 2));
std(:, 2) = sqrt(mean((Data(:, :, 2) - repmat(mu(:, 2), 1, n)).^2, 2));
std(:, 3) = sqrt(mean((Data(:, :, 3) - repmat(mu(:, 3), 1, n)).^2, 2));
std(:, 4) = sqrt(mean((Data(:, :, 4) - repmat(mu(:, 4), 1, n)).^2, 2));
std(:, 5) = sqrt(mean((Data(:, :, 5) - repmat(mu(:, 5), 1, n)).^2, 2));
std(:, 6) = sqrt(mean((Data(:, :, 6) - repmat(mu(:, 6), 1, n)).^2, 2));

clear hip_angle_l hip_angle_r knee_angle_l knee_angle_r ankle_angle_l ankle_angle_r

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
lambda = V(1:k);

eig_func = zeros(l, 6, k);
for j = 1:k
    eig_func(:, :, j) = reshape(U(:,j), [l, 6]);
end

score = zeros(n, k);
for i = 1:n
    sample = zeros(l, 6);
    sample(:, 1) = Data(:, i, 1);
    sample(:, 2) = Data(:, i, 2);
    sample(:, 3) = Data(:, i, 3);
    sample(:, 4) = Data(:, i, 4);
    sample(:, 5) = Data(:, i, 5);
    sample(:, 6) = Data(:, i, 6);
    sample = (sample - mu)./std;
    for j = 1:k
        score(i, j) = sum(sum(sample.*eig_func(:, :, j)));
    end
end

Data_re = zeros(l, 6, n);
for i = 1:n
    sample = zeros(l, 6);
    for j = 1:k
        sample = sample + score(i, j)*eig_func(:, :, j);
    end
    Data_re(:, :, i) = sample.*std + mu;
end

Data_re2 = zeros(l, 6, n);
for i = 1:n
    sample = zeros(l, 6);
    for j = 1:2
        sample = sample + score(i, j)*eig_func(:, :, j);
    end
    Data_re2(:, :, i) = sample.*std + mu;
end

%% plot 1

% yl = [-30 30];
% r = randi([1, n], 1, 60);
yl = [-140 120];
r = 1:n;
t = 0.01*((1:l)-1);

figure()
% clf

subplot(1, 3, 1)
hold on
for i = r
    plot(t, Data(:, i, 1), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data(:, i, 3), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data(:, i, 5), 'Color', '#EDB120', 'LineWidth', 0.8)
end
xlabel("time(sec)")
ylabel("angle(degree)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(1, 3, 2)
hold on
fill([t, flip(t)], [mu(:,1)' + 2*std(:,1)', flip(mu(:,1)' - 2*std(:,1)')], 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'LineStyle', 'none');
fill([t, flip(t)], [mu(:,3)' + 2*std(:,3)', flip(mu(:,3)' - 2*std(:,3)')], 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'LineStyle', 'none');
fill([t, flip(t)], [mu(:,5)' + 2*std(:,5)', flip(mu(:,5)' - 2*std(:,5)')], 1, 'FaceColor', '#EDB120', 'FaceAlpha', 0.2, 'LineStyle',' none');
plot(t, mu(:,1), 'Color', '#0072BD', 'LineWidth', 1.2);
plot(t, mu(:,3), 'Color', '#D95319', 'LineWidth', 1.2);
plot(t, mu(:,5), 'Color', '#EDB120', 'LineWidth', 1.2);
xlabel("time(sec)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(1, 3, 3)
hold on
for i = r
    plot(t, Data_re(:, 1, i), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data_re(:, 3, i), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data_re(:, 5, i), 'Color', '#EDB120', 'LineWidth', 0.8)
end
title(join(["( n_d = ", num2str(k), ', \eta = ', num2str(sum(lambda)/sum(V)), " )"]))
legend("Hip", "Knee", "Ankle")
xlabel("time(sec)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

% subplot(1,3,3)
% hold on
% plot([-500, 500], [0 0], 'k')
% plot([0, 0], [-500 500], 'k')
% plot(score(:,1), score(:,2), 'b.', 'MarkerSize', 6)
% title(join(["( n_d = ", num2str(k), ', \eta = ', num2str(sum(lambda)/sum(V)), " )"]))
% xlim(range)
% ylim(range)
% xlabel("PC1 Score")
% ylabel("PC2 Score")
% hold off
% box on
% grid on

%% plot 2

yl = [-120 80];
% range = [-30 30];
% r = randi([1, n], 1, 60);
r = 1:n;
t = 0.01*((1:l)-1);

figure()
% clf

subplot(2, 3, 1)
hold on
for i = r
    plot(t, Data(:, i, 1), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data(:, i, 3), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data(:, i, 5), 'Color', '#EDB120', 'LineWidth', 0.8)
end
ylabel("angle(degree)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(2, 3, 2)
hold on
fill([t, flip(t)], [mu(:,1)' + 2*std(:,1)', flip(mu(:,1)' - 2*std(:,1)')], 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'LineStyle', 'none');
fill([t, flip(t)], [mu(:,3)' + 2*std(:,3)', flip(mu(:,3)' - 2*std(:,3)')], 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'LineStyle', 'none');
fill([t, flip(t)], [mu(:,5)' + 2*std(:,5)', flip(mu(:,5)' - 2*std(:,5)')], 1, 'FaceColor', '#EDB120', 'FaceAlpha', 0.2, 'LineStyle',' none');
plot(t, mu(:,1), 'Color', '#0072BD', 'LineWidth', 1.2);
plot(t, mu(:,3), 'Color', '#D95319', 'LineWidth', 1.2);
plot(t, mu(:,5), 'Color', '#EDB120', 'LineWidth', 1.2);
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(2, 3, 3)
hold on
for i = r
    plot(t, Data_re(:, 1, i), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data_re(:, 3, i), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data_re(:, 5, i), 'Color', '#EDB120', 'LineWidth', 0.8)
end
title(join(["( n_d = ", num2str(k), ', \eta = ', num2str(sum(lambda)/sum(V)), " )"]))
legend("Hip", "Knee", "Ankle")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(2, 3, 4)
hold on
for i = r
    plot(t, Data(:, i, 2), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data(:, i, 4), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data(:, i, 6), 'Color', '#EDB120', 'LineWidth', 0.8)
end
xlabel("time(sec)")
ylabel("angle(degree)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(2, 3, 5)
hold on
fill([t, flip(t)], [mu(:,2)' + 2*std(:,2)', flip(mu(:,2)' - 2*std(:,2)')], 1, 'FaceColor', '#0072BD', 'FaceAlpha', 0.2, 'LineStyle', 'none');
fill([t, flip(t)], [mu(:,4)' + 2*std(:,4)', flip(mu(:,4)' - 2*std(:,4)')], 1, 'FaceColor', '#D95319', 'FaceAlpha', 0.2, 'LineStyle', 'none');
fill([t, flip(t)], [mu(:,6)' + 2*std(:,6)', flip(mu(:,6)' - 2*std(:,6)')], 1, 'FaceColor', '#EDB120', 'FaceAlpha', 0.2, 'LineStyle',' none');
plot(t, mu(:,2), 'Color', '#0072BD', 'LineWidth', 1.2);
plot(t, mu(:,4), 'Color', '#D95319', 'LineWidth', 1.2);
plot(t, mu(:,6), 'Color', '#EDB120', 'LineWidth', 1.2);
xlabel("time(sec)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

subplot(2, 3, 6)
hold on
for i = r
    plot(t, Data_re(:, 2, i), 'Color', '#0072BD', 'LineWidth', 0.8)
    plot(t, Data_re(:, 4, i), 'Color', '#D95319', 'LineWidth', 0.8)
    plot(t, Data_re(:, 6, i), 'Color', '#EDB120', 'LineWidth', 0.8)
end
xlabel("time(sec)")
xlim([t(1) t(end)])
ylim(yl)
hold off
box on
grid on

% subplot(2,3,6)
% hold on
% plot([-500, 500], [0 0], 'k')
% plot([0, 0], [-500 500], 'k')
% plot(score(:,1), score(:,2), 'b.', 'MarkerSize', 6)
% title(join(["( n_d = ", num2str(k), ', \eta = ', num2str(sum(lambda)/sum(V)), " )"]))
% xlim(range)
% ylim(range)
% xlabel("PC1 Score")
% ylabel("PC2 Score")
% hold off
% box on
% grid on

%% plot principal component

t = 0.01*((1:l)-1);

figure()
subplot(2, 3, 1)
hold on
plot(t, mu(:, 1), 'Color', '#0072BD', 'LineWidth', 1.2)
plot(t, mu(:, 3), 'Color', '#D95319', 'LineWidth', 1.2)
plot(t, mu(:, 5), 'Color', '#EDB120', 'LineWidth', 1.2)
title('Mean (Left Leg)')
xlim([t(1) t(end)])
ylim([-120, 80])
hold off
box on
grid on

subplot(2, 3, 2)
hold on
plot(t, mu(:, 2), 'Color', '#0072BD', 'LineWidth', 1.2)
plot(t, mu(:, 4), 'Color', '#D95319', 'LineWidth', 1.2)
plot(t, mu(:, 6), 'Color', '#EDB120', 'LineWidth', 1.2)
title('Mean (Right Leg)')
legend("Hip", "Knee", "Ankle")
xlim([t(1) t(end)])
ylim([-120, 80])
hold off
box on
grid on

subplot(2, 3, 4)
hold on
plot(t, std(:,1), 'Color', '#0072BD', 'LineWidth', 1.2);
plot(t, std(:,3), 'Color', '#D95319', 'LineWidth', 1.2);
plot(t, std(:,5), 'Color', '#EDB120', 'LineWidth', 1.2);
title('Std (Left Leg)')
xlim([t(1) t(end)])
ylim([0 10])
hold off
box on
grid on

subplot(2, 3, 5)
hold on
plot(t, std(:,2), 'Color', '#0072BD', 'LineWidth', 1.2);
plot(t, std(:,4), 'Color', '#D95319', 'LineWidth', 1.2);
plot(t, std(:,6), 'Color', '#EDB120', 'LineWidth', 1.2);
title('Std (Right Leg)')
xlim([t(1) t(end)])
ylim([0 10])
hold off
box on
grid on

figure()
for i = 1:6
    subplot(3, 4, 2*i - 1)
    hold on
    plot(t, eig_func(:, 1, i), 'Color', '#0072BD', 'LineWidth', 1.2)
    plot(t, eig_func(:, 3, i), 'Color', '#D95319', 'LineWidth', 1.2)
    plot(t, eig_func(:, 5, i), 'Color', '#EDB120', 'LineWidth', 1.2)
    title(join(["PC-", num2str(i), ', \lambda = ', num2str(lambda(i)), ', \eta = ', num2str(lambda(i)/sum(lambda))]), 'FontSize', 9)
    xlim([t(1) t(end)])
    ylim([-0.155 0.155])
    hold off
    box on
    grid on
    
    subplot(3, 4, 2*i)
    hold on
    plot(t, eig_func(:, 2, i), 'Color', '#0072BD', 'LineWidth', 1.2)
    plot(t, eig_func(:, 4, i), 'Color', '#D95319', 'LineWidth', 1.2)
    plot(t, eig_func(:, 6, i), 'Color', '#EDB120', 'LineWidth', 1.2)
    title(join(["PC-", num2str(i), ', \lambda = ', num2str(lambda(i)), ', \eta = ', num2str(lambda(i)/sum(lambda))]), 'FontSize', 9)
    xlim([t(1) t(end)])
    ylim([-0.155 0.155])
    hold off
    box on
    grid on
end

figure()
subplot(3, 4, 1)
hold on
plot(t, eig_func(:, 1, 7), 'Color', '#0072BD', 'LineWidth', 1.2)
plot(t, eig_func(:, 3, 7), 'Color', '#D95319', 'LineWidth', 1.2)
plot(t, eig_func(:, 5, 7), 'Color', '#EDB120', 'LineWidth', 1.2)
title(join(["PC-", num2str(7), ', \lambda = ', num2str(lambda(7)), ', \eta = ', num2str(lambda(7)/sum(lambda))]), 'FontSize', 9)
xlim([t(1) t(end)])
ylim([-0.155 0.155])
hold off
box on
grid on

subplot(3, 4, 2)
hold on
plot(t, eig_func(:, 2, 7), 'Color', '#0072BD', 'LineWidth', 1.2)
plot(t, eig_func(:, 4, 7), 'Color', '#D95319', 'LineWidth', 1.2)
plot(t, eig_func(:, 6, 7), 'Color', '#EDB120', 'LineWidth', 1.2)
title(join(["PC-", num2str(7), ', \lambda = ', num2str(lambda(7)), ', \eta = ', num2str(lambda(7)/sum(lambda))]), 'FontSize', 9)
xlim([t(1) t(end)])
ylim([-0.155 0.155])
hold off
box on
grid on

%%

figure()
plot(score(:, 1), score(:, 2), '.', 'MarkerSize', 10)
axis([-80, 80, -80, 80])
box on
grid on


