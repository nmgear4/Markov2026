N = 1e4;

samples = zeros(N,1);
n = 0;
proposals = 0;

lambda = 0.2;
c = 1/(0.16*exp(1));

tic

while n < N
        U1 = rand;
        X = -log(U1)/lambda;
        U2 = rand;
        proposals = proposals + 1;

        if U2 < X*exp(-X)/(c*lambda*exp(-lambda*X))
            n = n+1;
            samples(n) = X;
        end
end


elapsed_t = toc;

acceptance_frac = N/proposals;
mean_t = elapsed_t/N;

fprintf('Acceptance fraction = %.4f\n', acceptance_frac);
fprintf('Theoretical acceptance = %.4f\n', exp(1)/4);
fprintf('Mean time per accepted sample = %.8f seconds\n', mean_t);

% Histogram + true destiny
histogram(samples,'Normalization','pdf');
hold on

x = linspace(0,15,500);
plot(x,x.*exp(-x),'r','LineWidth',2);

xlabel('x');
ylabel('Density');
legend('Simulated samples','f(x) = xe^{-x}');
title('\lambda = 0.5');
grid on
