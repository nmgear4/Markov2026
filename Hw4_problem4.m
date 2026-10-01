P = [0   1/2 1/2 0   0   0;
     1/4 0   0   1/2 0   1/4;
     3/4 0   0   0   0   1/4;
     0   0   0   3/4 1/4 0;
     0   0   0   1/2 1/2 0;
     0   0   0   0   0   1];

Nsim = 10^4;

state_names = {'U','I','M'};

h_hat = zeros(3,1);
g_hat = zeros(3,1);
tauF_hat = zeros(3,1);
tauA_hat = zeros(3,1);

all_T = cell(3,1);
all_fate = cell(3,1);

for start = 1:3

    T = zeros(Nsim,1);
    fate = zeros(Nsim,1);

    for sim = 1:Nsim

        state = start;
        t = 0;

        while state <= 3

            t = t + 1;

            r = rand;

            cumulative = cumsum(P(state,:));

            next_state = find(r <= cumulative,1);

            state = next_state;

        end

        T(sim) = t;

        if state == 4 || state == 5
            fate(sim) = 1;
        else
            fate(sim) = 2;
        end

    end

    all_T{start} = T;
    all_fate{start} = fate;

    h_hat(start) = mean(fate == 1);

    g_hat(start) = mean(T);

    T_F = T(fate == 1);
    tauF_hat(start) = mean(T_F);

    T_A = T(fate == 2);
    tauA_hat(start) = mean(T_A);

end

h_exact = [1/2; 5/8; 3/8];
g_exact = [4; 2; 4];
tauF_exact = [4; 9/5; 5];
tauA_exact = [4; 7/3; 17/5];

fprintf('\n');
fprintf('Simulation results, N = %d\n',Nsim);
fprintf('Start     h_hat      h_exact    g_hat      g_exact\n');

for x = 1:3

    fprintf('%s       %.5f     %.5f     %.5f     %.5f\n', ...
        state_names{x}, ...
        h_hat(x), ...
        h_exact(x), ...
        g_hat(x), ...
        g_exact(x));

end

fprintf('\n');
fprintf('Conditional mean first passage times\n');
fprintf('Start     tauF_hat   tauF_exact   tauA_hat   tauA_exact\n');

for x = 1:3

    fprintf('%s       %.5f      %.5f       %.5f      %.5f\n', ...
        state_names{x}, ...
        tauF_hat(x), ...
        tauF_exact(x), ...
        tauA_hat(x), ...
        tauA_exact(x));

end



T_I = all_T{2};
fate_I = all_fate{2};

T_fold = T_I(fate_I == 1);
T_aggregate = T_I(fate_I == 2);

Q = P(1:3,1:3);

R = [P(1,4)+P(1,5), P(1,6);
    P(2,4)+P(2,5), P(2,6);
    P(3,4)+P(3,5), P(3,6)];

maxT = max(T_I);

n_values = 1:maxT;

pmf_F = zeros(size(n_values));
pmf_A = zeros(size(n_values));

for k = 1:length(n_values)

    n = n_values(k);

    if n == 1
        Qpower = eye(3);
    else
        Qpower = Q^(n-1);
    end

    B = Qpower * R;

    prob_F = B(2,1);
    prob_A = B(2,2);

    pmf_F(k) = prob_F / h_exact(2);
    pmf_A(k) = prob_A / (1-h_exact(2));

end

figure;

histogram(T_fold, ...
    'Normalization','probability', ...
    'BinMethod','integers', ...
    'FaceColor',[0.2 0.5 0.9], ...
    'EdgeColor','none');

hold on;

stem(n_values,pmf_F, ...
    'r','LineWidth',1.5, ...
    'Marker','none');

xlabel('T');
ylabel('Probability');
title('Starting from I: Folded Runs');

legend('Simulation','Exact conditional PMF', ...
    'Location','best');

grid on;

figure;

histogram(T_aggregate, ...
    'Normalization','probability', ...
    'BinMethod','integers', ...
    'FaceColor',[0.9 0.5 0.2], ...
    'EdgeColor','none');

hold on;

stem(n_values,pmf_A, ...
    'b','LineWidth',1.5, ...
    'Marker','none');

xlabel('T');
ylabel('Probability');
title('Starting from I: Aggregated Runs');

legend('Simulation','Exact conditional PMF', ...
    'Location','best');

grid on;

figure;

histogram(T_fold, ...
    'Normalization','probability', ...
    'BinMethod','integers', ...
    'FaceColor',[0.2 0.5 0.9], ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none');

hold on;

histogram(T_aggregate, ...
    'Normalization','probability', ...
    'BinMethod','integers', ...
    'FaceColor',[0.9 0.5 0.2], ...
    'FaceAlpha',0.5, ...
    'EdgeColor','none');

stem(n_values,pmf_F, ...
    'b','LineWidth',1.5, ...
    'Marker','none');

stem(n_values,pmf_A, ...
    'r','LineWidth',1.5, ...
    'Marker','none');

xlabel('T');
ylabel('Probability');
title('Starting from I: Conditional Distribution of T');

legend('Folded simulation', ...
    'Aggregated simulation', ...
    'Exact folded PMF', ...
    'Exact aggregated PMF', ...
    'Location','best');

grid on;
