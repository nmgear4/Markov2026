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

