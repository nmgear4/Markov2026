lambda = 4;
H = 2e5;

Nexp = zeros(H,1);
firstExp = nan(H,1);
halfExp = zeros(H,1);

for j = 1:H
    time = 0;
    count = 0;
    first = NaN;
    countFirstHalf = 0;

    while true
        time = time - log(rand)/lambda;
        if time > 1
            break;
        end

        count = count + 1;

        if count == 1
            first = time;
        end

        if time <= 0.5
            countFirstHalf = countFirstHalf + 1;
        end
    end

    Nexp(j) = count;
    firstExp(j) = first;
    halfExp(j) = countFirstHalf;
end

Nuni = poissrnd(lambda,H,1);
firstUni = nan(H,1);
halfUni = zeros(H,1);

for j = 1:H
    K = Nuni(j);
    if K > 0
        points = sort(rand(K,1));
        firstUni(j) = points(1);
        halfUni(j) = sum(points <= 0.5);
    end
end

k = 0:14;
poisPMF = exp(-lambda + k*log(lambda) ...
    - gammaln(k+1));

fprintf('One-hour count summaries:\n');
fprintf('Method       Mean       Variance\n');
fprintf('Gaps       %8.5f    %8.5f\n', ...
    mean(Nexp),var(Nexp,1));
fprintf('Uniform    %8.5f    %8.5f\n', ...
    mean(Nuni),var(Nuni,1));

figure('Color','w');
tiledlayout(1,2);

nexttile;
histogram(Nexp,'BinMethod','integers', ...
    'Normalization','probability');
hold on;
plot(k,poisPMF,'ro-','LineWidth',1.4);
xlim([-0.5 12.5]);
xlabel('N(1)');
ylabel('Probability');
title('Exponential gaps');
legend('Simulation','Poisson(4) pmf');
grid on;
hold off;

nexttile;
histogram(Nuni,'BinMethod','integers', ...
    'Normalization','probability');
hold on;
plot(k,poisPMF,'ro-','LineWidth',1.4);
xlim([-0.5 12.5]);
xlabel('N(1)');
ylabel('Probability');
title('Poisson count + uniform points');
legend('Simulation','Poisson(4) pmf');
grid on;
hold off;

R = 10;
duration = 10;

figure('Color','w');
tiledlayout(2,1);

nexttile;
hold on;
for r = 1:R
    times = rasterGaps(lambda,duration);
    scatter(times,r*ones(size(times)),12,'filled');
end
xlim([0 duration]);
ylim([0 R+1]);
ylabel('Realization');
title('Raster: exponential gaps');
grid on;
hold off;

nexttile;
hold on;
for r = 1:R
    times = rasterUniform(lambda,duration);
    scatter(times,r*ones(size(times)),12,'filled');
end
xlim([0 duration]);
ylim([0 R+1]);
xlabel('Time (hours)');
ylabel('Realization');
title('Raster: Poisson counts + uniform points');
grid on;
hold off;

T3 = sum(-log(rand(H,3))/lambda,2);

t = linspace(0,3,600);
gammaDensity = lambda^3 .* t.^2 .* exp(-lambda*t)/factorial(2);

figure('Color','w');
histogram(T3,100,'Normalization','pdf', ...
    'FaceColor',[0.65 0.78 0.92], ...
    'EdgeColor','none');
hold on;
plot(t,gammaDensity,'r-','LineWidth',2);
xlabel('T_3 (hours)');
ylabel('Density');
title('Third arrival time: Gamma(3,4)');
legend('Simulation','Exact Gamma density');
grid on;
hold off;

fprintf('\nT3 simulation: mean = %.5f; theory = %.5f\n', ...
    mean(T3),3/lambda);
fprintf('T3 simulation: variance = %.5f; theory = %.5f\n', ...
    var(T3,1),3/lambda^2);

idxE = Nexp == 3;
idxU = Nuni == 3;

pHalfZeroE = mean(halfExp(idxE)==0);
pHalfZeroU = mean(halfUni(idxU)==0);

meanFirstE = mean(firstExp(idxE));
meanFirstU = mean(firstUni(idxU));

fprintf('\nConditional on N(1)=3:\n');
fprintf('Statistic                 Gaps       Uniform\n');
fprintf('P(N(1/2)=0 | N(1)=3)   %8.5f   %8.5f\n', ...
    pHalfZeroE,pHalfZeroU);
fprintf('E[T1 | N(1)=3]         %8.5f   %8.5f\n', ...
    meanFirstE,meanFirstU);
fprintf('Theory                  %8.5f   %8.5f\n', ...
    1/8,1/4);

secondExp = Nexp-halfExp;
secondUni = Nuni-halfUni;

Cexp = corrcoef(halfExp,secondExp);
Cuni = corrcoef(halfUni,secondUni);

fprintf('\nCorrelation between half-hour counts:\n');
fprintf('Gaps method:   %.6f\n',Cexp(1,2));
fprintf('Uniform method: %.6f\n',Cuni(1,2));
fprintf('Theory:         0\n');

figure('Color','w');
histogram(firstExp(idxE),30,'Normalization','pdf');
hold on;
histogram(firstUni(idxU),30,'Normalization','pdf');
s = linspace(0,1,300);
plot(s,3*(1-s).^2,'k-','LineWidth',2);
xlabel('First arrival time (hours)');
ylabel('Density');
title('First arrival, conditional on N(1)=3');
legend('Gaps','Uniform points','Conditional density');
grid on;
hold off;

mValues = 10:10:2000;
dt = 1./mValues;
TV = zeros(size(mValues));

for j = 1:numel(mValues)
    m = mValues(j);
    p = lambda/m;

    k = 0:m;

    logBin = gammaln(m+1) ...
        - gammaln(k+1) ...
        - gammaln(m-k+1) ...
        + k*log(p) ...
        + (m-k)*log1p(-p);

    binPMF = exp(logBin);

    poisPMF = exp(-lambda + k*log(lambda) ...
        - gammaln(k+1));

    poissonTail = max(0,1-sum(poisPMF));

    TV(j) = 0.5*(sum(abs(binPMF-poisPMF)) + poissonTail);
end

fitCoeffs = polyfit(log(dt),log(TV),1);
slope = fitCoeffs(1);
intercept = fitCoeffs(2);

fprintf('Fitted log-log slope = %.5f\n',slope);
fprintf('Theoretical asymptotic slope is 1.\n');

figure('Color','w');
loglog(dt,TV,'b.','MarkerSize',12);
hold on;
loglog(dt,exp(intercept)*dt.^slope,'r-','LineWidth',2);
set(gca,'XDir','reverse');
xlabel('\Delta t');
ylabel('Total variation distance');
title(sprintf('Binomial approximation; fitted slope = %.3f',slope));
legend('Computed TV distance','Log-log regression', ...
    'Location','best');
grid on;
hold off;

function times = rasterGaps(lambda,duration)
    times = [];
    time = 0;
    while true
        time = time - log(rand)/lambda;
        if time > duration
            break;
        end
        times(end+1) = time;
    end
end

function times = rasterUniform(lambda,duration)
    numHours = ceil(duration);
    times = [];

    for h = 0:numHours-1
        K = poissrnd(lambda);
        points = h + rand(K,1);
        times = [times; points(points < duration)];
    end

    times = sort(times(:)).';
end
