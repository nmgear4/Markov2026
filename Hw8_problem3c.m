T = 720;
lambda_max = 20;

Nproposal = poissrnd(lambda_max*T);

U = T*rand(Nproposal,1);

lambda_U = 12 + 8*cos(pi*(U-16)/12);

keep = rand(Nproposal,1) < lambda_U/lambda_max;

arrival = sort(U(keep));

fraction_kept = length(arrival)/Nproposal;

fprintf('Number of proposals: %d\n', Nproposal);
fprintf('Number retained: %d\n', length(arrival));
fprintf('Fraction retained: %.4f\n', fraction_kept);

hourly_counts = zeros(30,24);

for d = 1:30
    for h = 0:23

        lower = 24*(d-1) + h;
        upper = lower + 1;

        hourly_counts(d,h+1) = ...
            sum(arrival >= lower & arrival < upper);

    end
end

mean_hourly = mean(hourly_counts,1);

true_hourly = zeros(1,24);

for h = 0:23

    true_hourly(h+1) = ...
        12 + (96/pi)*( ...
        sin(pi*(h+1-16)/12) ...
        - sin(pi*(h-16)/12));

end

figure;

plot(0:23, mean_hourly, 'bo-', ...
    'LineWidth', 1.5, 'MarkerSize', 6);

hold on;

plot(0:23, true_hourly, 'r-', ...
    'LineWidth', 2);

xlabel('Hour after midnight');
ylabel('Mean number of arrivals per hour');

title('Average hourly arrivals over 30 days');

legend('Simulated mean', 'True expected count', ...
    'Location', 'best');

xticks(0:23);
grid on;

early_morning = zeros(30,1);

for d = 1:30

    lower = 24*(d-1) + 3;
    upper = lower + 1;

    early_morning(d) = ...
        sum(arrival >= lower & arrival < upper);

end

sample_mean = mean(early_morning);
sample_variance = var(early_morning);

mu = 12 - (96/pi)*sin(pi/12);

fprintf('\n3 AM TO 4 AM RESULTS\n');

fprintf('Simulated mean: %.4f\n', sample_mean);
fprintf('Simulated variance: %.4f\n', sample_variance);

fprintf('Theoretical mean: %.4f\n', mu);
fprintf('Theoretical variance: %.4f\n', mu);

n = length(arrival);

negloglik = @(theta) objective(theta, arrival, T);

initial = [10, 5, 14];

options = optimset('Display','iter', ...
    'MaxIter',5000, ...
    'MaxFunEvals',10000, ...
    'TolX',1e-9, ...
    'TolFun',1e-9);

theta_hat = fminsearch(negloglik, initial, options);

a_hat = theta_hat(1);
b_hat = theta_hat(2);
phi_hat = theta_hat(3);

LL_NHPP = -negloglik(theta_hat);

fprintf('\nMAXIMUM LIKELIHOOD ESTIMATES\n');

fprintf('Estimated a = %.5f\n', a_hat);
fprintf('Estimated b = %.5f\n', b_hat);
fprintf('Estimated phi = %.5f\n', phi_hat);

fprintf('Maximum log-likelihood = %.5f\n', LL_NHPP);

lambda_constant = n/T;

LL_constant = n*log(lambda_constant) ...
    - lambda_constant*T;

fprintf('\nCONSTANT RATE MODEL\n');

fprintf('Estimated constant rate = %.5f\n', ...
    lambda_constant);

fprintf('Constant model log-likelihood = %.5f\n', ...
    LL_constant);

fprintf('\nMODEL COMPARISON\n');

fprintf('NHPP log-likelihood = %.5f\n', LL_NHPP);
fprintf('Constant log-likelihood = %.5f\n', LL_constant);


function value = objective(theta, arrival, T)

a = theta(1);
b = theta(2);
phi = theta(3);

if a <= abs(b) || a <= 0
    value = 1e10;
    return;
end

lambda = a + b*cos(pi*(arrival-phi)/12);

if any(lambda <= 0)
    value = 1e10;
    return;
end

Lambda = a*T + (12*b/pi)*( ...
    sin(pi*(T-phi)/12) ...
    - sin(-pi*phi/12));

value = -sum(log(lambda)) + Lambda;

end
