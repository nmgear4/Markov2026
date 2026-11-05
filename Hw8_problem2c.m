N = 2e5;
H = 5;
q = 0.4;
lambda_s = 3;
lambda_b = 2;

mu_s = q * lambda_s;
mu_b = q * lambda_b;

star_original = poissrnd(lambda_s * H, N, 1);
sky_original  = poissrnd(lambda_b * H, N, 1);

star1 = binornd(star_original, q);
sky1  = binornd(sky_original, q);

total1 = star1 + sky1;


star2 = poissrnd(mu_s * H, N, 1);
sky2  = poissrnd(mu_b * H, N, 1);

total2 = star2 + sky2;


fprintf('\nMETHOD 1: THINNING\n');

fprintf('Mean star counts: %.4f\n', mean(star1));
fprintf('Variance star counts: %.4f\n', var(star1));

fprintf('Mean sky counts: %.4f\n', mean(sky1));
fprintf('Variance sky counts: %.4f\n', var(sky1));

fprintf('Correlation (star, sky): %.4f\n', ...
    corr(star1, sky1));

fprintf('Mean total counts: %.4f\n', mean(total1));
fprintf('Variance total counts: %.4f\n', var(total1));

fprintf('\nMETHOD 2: DIRECT GENERATION\n');

fprintf('Mean star counts: %.4f\n', mean(star2));
fprintf('Variance star counts: %.4f\n', var(star2));

fprintf('Mean sky counts: %.4f\n', mean(sky2));
fprintf('Variance sky counts: %.4f\n', var(sky2));

fprintf('Correlation (star, sky): %.4f\n', ...
    corr(star2, sky2));

fprintf('Mean total counts: %.4f\n', mean(total2));
fprintf('Variance total counts: %.4f\n', var(total2));

fprintf('\nTOTAL COUNT COMPARISON\n');

fprintf('Theoretical Poisson mean: 10\n');
fprintf('Theoretical Poisson variance: 10\n');

fprintf('Method 1 mean: %.4f\n', mean(total1));
fprintf('Method 1 variance: %.4f\n', var(total1));

fprintf('Method 2 mean: %.4f\n', mean(total2));
fprintf('Method 2 variance: %.4f\n', var(total2));

idx1 = (total1 == 10);
idx2 = (total2 == 10);

prob1 = mean(star1(idx1) == 6);
prob2 = mean(star2(idx2) == 6);

fprintf('\nCONDITIONAL PROBABILITY\n');

fprintf('Method 1 P(6 star | 10 total): %.5f\n', prob1);
fprintf('Method 2 P(6 star | 10 total): %.5f\n', prob2);

fprintf('Theoretical probability: %.5f\n', ...
    nchoosek(10,6)*0.6^6*0.4^4);

Ts = -log(rand(N,1)) / mu_s;
Tb = -log(rand(N,1)) / mu_b;

T = min(Ts,Tb);

M = max(Ts,Tb);

starFirst = (Ts < Tb);


P_star_first = mean(starFirst);

P_star_first_1s = mean(starFirst & (T <= 1));

mean_M = mean(M);

correlation_TM = corr(T, double(starFirst));

fprintf('\nFIRST DETECTION TIME RESULTS\n');

fprintf('P(star first): %.5f\n', P_star_first);

fprintf('P(star first and within 1 second): %.5f\n', ...
    P_star_first_1s);

fprintf('Expected time until both types appear: %.5f\n', ...
    mean_M);

fprintf('Correlation between first arrival time and type: %.5f\n', ...
    correlation_TM);
