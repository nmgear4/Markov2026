N = 1e5;

U1 = rand(N,1);
U2 = rand(N,1);

T = zeros(N,1);

fast = (U1 < 0.9);

T(fast) = -log(U2(fast))/1000;
T(~fast) = -log(U2(~fast))/10;

binWidth = 0.001;

histogram(T, 'BinWidth', binWidth, 'Normalization', 'pdf');
hold on

t = linspace(0,1,1000);
f = 0.9*1000*exp(-1000*t) + 0.1*10*exp(-10*t);

empirical_mean = mean(T);

fprintf('Empirical mean = %.6f s (%.3f ms)\n',empirical_mean, 1000*empirical_mean);

empirical_tail = mean(T > 0.05);
theoretical_tail = 0.9*exp(-1000*0.05) + 0.1*exp(-10*0.05);

fprintf('Empirical P(T > 50 ms) = %.5f\n', empirical_tail);
fprintf('Theoretical P(T > 50 ms) = %.5f\n', theoretical_tail);

plot(t, f, 'r', 'LineWidth', 2);
    set(gca, 'YScale', 'log');
    xlabel('Dwell time t(s)');
    ylabel('Density');
    title('Mixture of Exponential Dwell Times');
    legend('Empirical PDF', 'Theoretical PDF');
    grid on
