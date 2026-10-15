N = 100;
K = 10;
gamma = 1;
dt = 0.005;

R0vals = [0.5, 0.8, 1.5, 3];

fprintf('Expected absorption times\n');
fprintf(' R0        e1 (steps)       e1*dt (time units)\n');

for R0 = R0vals
    beta = R0*gamma/N;

    n = (1:K-1)';
    p = beta*n.*(N-n)*dt;
    q = gamma*n*dt;

    M = diag(p+q);

    M = M + diag(-p(1:end-1),1);

    M = M + diag(-q(2:end),-1);

    e = M\ones(K-1,1);

    fprintf('%4.1f     %12.4f       %12.6f\n',...
        R0,e(1),e(1)*dt);
end

R0 = 1.5;
beta = R0*gamma/N;

numTrials = 20000;
totalSteps = zeros(numTrials,1);

rng(1);

for trial = 1:numTrials
    n = 1;
    steps = 0;

    while n > 0 && n < K
        p = beta*n*(N-n)*dt;
        q = gamma*n*dt;
        u = rand;

        if u < p
            n = n + 1;
        elseif u < p+q
            n = n - 1;
        end

        steps = steps + 1;
    end

    totalSteps(trial) = steps;
end

simTime = mean(totalSteps)*dt;
seTime = std(totalSteps)/sqrt(numTrials)*dt;

fprintf('\nSimulation check, R0 = %.1f\n',R0);
fprintf('Estimated E[T]*dt = %.6f\n',simTime);
fprintf('Monte Carlo standard error = %.6f\n',seTime);
