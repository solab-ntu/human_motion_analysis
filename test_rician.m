%% rician distribution

n = 3000000;
x1 = normrnd(1, 5, [1 n]);
x2 = normrnd(1, 5, [1 n]);
d = sqrt(x1.^2 + x2.^2);
d = sort(d);

std = rms(d);
display([std, 2*std, sum(d <= 2*std)/length(d)])

%%

figure()
subplot(1, 2, 1)
hold on
plot(x1, x2, '.', 'markerSize', 5)
plot(5, 5, '.', 'markerSize', 12)
plot([-10 20], [0 0], 'k-')
plot([0 0], [-10 20], 'k-')
axis([-10 20 -10 20])
hold off
box on

ub = 30;
pd = fitdist(d', 'rician');
xx = 0:0.5:ub;
yy = pdf(pd, xx);
subplot(1, 2, 2)
hold on
h = histogram(d, 'BinWidth' , 0.5, 'BinLimits', [0 ub]);
values = h.Values;
opt = optimoptions('fmincon','Display','off');
a = fmincon(@(a) sum((a*yy(2:end) - values).^2), 1, [], [], [], [], 0, Inf, [], opt);
plot(xx, a*yy, 'LineWidth', 1.0)
xlim([0 ub])
hold off
box on