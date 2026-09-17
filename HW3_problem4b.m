clear;
clc;
rng(1);

d0 = 10;
R = 2e4;
T = 1e4;

t = (1:T)';

x1 = zeros(R,1);
y1 = d0 * ones(R,1);
alive1 = true(R,1);

S1 = zeros(T,1);

for n = 1:T
    lamb_step = 2*(rand(R,1) < 0.5) - 1;
    lion_step = 2*(rand(R,1) < 0.5) - 1;

    x1(alive1) = x1(alive1) + lamb_step(alive1);
    y1(alive1) = y1(alive1) + lion_step(alive1);

    captured = alive1 & (x1 == y1);
    alive1(captured) = false;

    S1(n) = mean(alive1);
end

S1_cont = erf(d0 ./ (2*sqrt(t)));

fit_idx = (t >= 1e2 & t <= 1e4);
p1 = polyfit(log(t(fit_idx)), log(S1(fit_idx)), 1);
beta1 = -p1(1);

fprintf('Estimated beta1 = %.4f over [10^2, 10^4]\n', beta1);

figure; hold on;
loglog(t, S1, 'b', 'LineWidth', 1.5);
loglog(t, S1_cont, 'r--', 'LineWidth', 1.5);
legend('S1(t) simulation', 'S1(t) continuum');
xlabel('t'); ylabel('S1(t)');
title(sprintf('N=1 survival, beta1 = %.4f', beta1));
grid on;
