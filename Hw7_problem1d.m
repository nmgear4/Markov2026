alphaTrue = 1.5;
betaTrue = 0.8;
sampleSizes = [10 100 1000];
numSets = 2000;

scatterMean = zeros(size(sampleSizes));
scatterMedian = zeros(size(sampleSizes));

allMeans = cell(size(sampleSizes));
allMedians = cell(size(sampleSizes));

for j = 1:numel(sampleSizes)
    n = sampleSizes(j);

    U = rand(numSets,n);
    X = alphaTrue + betaTrue*tan(pi*(U - 0.5));

    estimatesMean = mean(X,2);
    estimatesMedian = median(X,2);

    allMeans{j} = estimatesMean;
    allMedians{j} = estimatesMedian;

    qMean = quantile(estimatesMean,[0.25 0.75]);
    qMedian = quantile(estimatesMedian,[0.25 0.75]);

    scatterMean(j) = (qMean(2)-qMean(1))/1.349;
    scatterMedian(j) = (qMedian(2)-qMedian(1))/1.349;

    fprintf(['n = %4d: mean scatter = %.5f, ' ...
        'median scatter = %.5f\n'], ...
        n,scatterMean(j),scatterMedian(j));
end

Results = table(sampleSizes(:),scatterMean(:),scatterMedian(:), ...
    'VariableNames',{'n','MeanScatter','MedianScatter'});
disp(Results);

figure('Color','w');
loglog(sampleSizes,scatterMean,'o-','LineWidth',1.5, ...
    'MarkerSize',7);
hold on;
loglog(sampleSizes,scatterMedian,'s-','LineWidth',1.5, ...
    'MarkerSize',7);
xlabel('Sample size n');
ylabel('IQR / 1.349');
title('Sampling variability of the mean and median');
legend('Sample mean','Sample median','Location','best');
grid on;
hold off;

figure('Color','w');
histogram(allMeans{3},60,'Normalization','pdf');
hold on;
histogram(allMedians{3},40,'Normalization','pdf');
xline(alphaTrue,'k--','True \alpha','LineWidth',1.5);
xlabel('Estimator');
ylabel('Density');
title('Sampling distributions for n = 1000');
legend('Sample mean','Sample median','True location');
grid on;
hold off;
