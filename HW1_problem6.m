% HW1 Problem 6b

clear; clc; close all;

% Parameters
Nvals = round(logspace(2,7,29));
Ptrue = 0.3763;

Nmax = max(Nvals);

rng(2026);

X = rand(Nmax,1);
Y = rand(Nmax,1);
Z = rand(Nmax,1);

success = (X.^2 + Y.^2 < Z) & (Z.^2 > X.*Y);

Phat = zeros(size(Nvals));

for k = 1:length(Nvals)
    N = Nvals(k);
    Phat(k) = mean(success(1:N));
end

%Plot
figure;
semilogx(Nvals, Phat, 'o-', 'LineWidth', 1.5);
hold on;
yline(Ptrue, 'r--', 'LineWidth', 1.5);

xlabel('Sample size N');
ylabel('Monte Carlo porbability estimate');
title('Monte Carlo estimate HW 1 6b');
legend('Monte Carlo estimate', 'Analytic value', ...
    'Location', "best");
grid on;
xlim([min(Nvals), max(Nvals)]);
