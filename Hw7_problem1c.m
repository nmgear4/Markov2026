x = [-1.77 2.47 1.34 2.18 13.05 1.60 1.50 ...
     -1.97 0.79 1.50 2.00 2.63 1.19 -2.31 0.87];
n = numel(x);
T = 5e4;
steps = [0.02 0.5 5];
nsteps = numel(steps);

logg = @(a,b) n*b - sum(log((exp(b).^2) + (x-a).^2));

chains = cell(nsteps,1);
acceptance = zeros(nsteps,1);
iactAlpha = zeros(nsteps,1);

for j = 1:nsteps
    s = steps(j);
    z = [0; 0];
    chain = zeros(T,2);
    accepted = 0;

    for t = 1:T
        y = z + s*randn(2,1);

        logRatio = logg(y(1),y(2)) - logg(z(1),z(2));

        if log(rand) < min(0,logRatio)
            z = y;
            accepted = accepted + 1;
        end

        chain(t,:) = z.';
    end

    chains{j} = chain;
    acceptance(j) = accepted/T;
    iactAlpha(j) = integratedTime(chain(:,1));

    fprintf(['s = %.2f: acceptance = %.4f, ' ...
             'IACT(alpha) = %.2f\n'], ...
             s, acceptance(j), iactAlpha(j));
end

figure('Color','w');
tiledlayout(3,1);

for j = 1:nsteps
    nexttile;
    plot(1:400,chains{j}(1:400,1),'b-');
    ylabel('\alpha');
    title(sprintf('s = %.2f',steps(j)));
    grid on;
    if j == 3
        xlabel('Iteration');
    end
end

[~,best] = min(iactAlpha);
sBest = steps(best);
chain = chains{best};

fprintf('\nBest step size by estimated IACT: %.2f\n',sBest);

post = chain(5001:end,:);
alphaChain = post(:,1);
betaChain = exp(post(:,2));

chainMean = [mean(alphaChain), mean(betaChain)];
chainCI = [quantile(alphaChain,[0.05 0.95]); ...
           quantile(betaChain,[0.05 0.95])];

fprintf('\nPosterior summaries from MCMC:\n');
fprintf('alpha mean = %.5f, 90%% interval = [%.5f, %.5f]\n', ...
    chainMean(1),chainCI(1,1),chainCI(1,2));
fprintf('beta  mean = %.5f, 90%% interval = [%.5f, %.5f]\n', ...
    chainMean(2),chainCI(2,1),chainCI(2,2));

aGrid = linspace(-5,7,401);
bGrid = linspace(log(0.05),log(15),401);
[A,B] = meshgrid(aGrid,bGrid);

logG = n*B;
for i = 1:n
    logG = logG - log(exp(2*B) + (x(i)-A).^2);
end

W = exp(logG - max(logG(:)));
W = W / sum(W(:));

betaGrid = exp(B);
meanAlphaGrid = sum(A(:).*W(:));
meanBetaGrid = sum(betaGrid(:).*W(:));

margA = sum(W,1);
margB = sum(W,2).';

cdfA = cumsum(margA);
cdfB = cumsum(margB);

qa = [interp1(cdfA,aGrid,0.05), ...
      interp1(cdfA,aGrid,0.95)];
qb = [interp1(cdfB,bGrid,0.05), ...
      interp1(cdfB,bGrid,0.95)];

qb = exp(qb);

fprintf('\nPosterior summaries from grid:\n');
fprintf('alpha mean = %.5f, 90%% interval = [%.5f, %.5f]\n', ...
    meanAlphaGrid,qa(1),qa(2));
fprintf('beta  mean = %.5f, 90%% interval = [%.5f, %.5f]\n', ...
    meanBetaGrid,qb(1),qb(2));

figure('Color','w');
contour(A,B,W,15,'LineWidth',1.1);
hold on;
plot(post(:,1),post(:,2),'.', ...
     'Color',[0.85 0.2 0.2],'MarkerSize',3);
xlabel('\alpha (km)');
ylabel('b = log(\beta)');
title(sprintf('Posterior contours and MCMC samples, s = %.2f',sBest));
legend('Grid posterior','MCMC samples','Location','best');
grid on;
hold off;

fprintf('\nData summaries:\n');
fprintf('Sample mean   = %.5f km\n',mean(x));
fprintf('Sample median = %.5f km\n',median(x));

function tau = integratedTime(v)
    v = v(:);
    v = v - mean(v);
    N = length(v);

    nfft = 2^nextpow2(2*N-1);
    F = fft(v,nfft);
    acov = real(ifft(F.*conj(F)));
    acov = acov(1:N);
    acov = acov/acov(1);

    tau = 1;
    k = 2;
    while k < N
        pairSum = acov(k) + acov(k+1);
        if pairSum <= 0
            break;
        end
        tau = tau + 2*pairSum;
        k = k + 2;
    end
end
