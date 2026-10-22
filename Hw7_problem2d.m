lambda = [3 2 1.5 1.2 1 0.75];
n = numel(lambda);
N = 2e5;
numGroups = 2^n;

E = zeros(numGroups,1);

PLast = zeros(numGroups,n);

for mask = 1:numGroups-1
    members = find(bitget(mask,1:n));
    L = sum(lambda(members));

    if numel(members) == 1
        i = members(1);
        E(mask+1) = 1/lambda(i);
        PLast(mask+1,i) = 1;
    else
        E(mask+1) = 1/L;
        for i = members
            smallerMask = bitset(mask,i,0);
            probNext = lambda(i)/L;

            E(mask+1) = E(mask+1) + ...
                probNext*E(smallerMask+1);

            PLast(mask+1,:) = PLast(mask+1,:) + ...
                probNext*PLast(smallerMask+1,:);
        end
    end
end

fullMask = numGroups-1;
Erecursive = E(fullMask+1);
pLastExact = PLast(fullMask+1,:);

Ealternating = 0;

for mask = 1:numGroups-1
    members = find(bitget(mask,1:n));
    L = sum(lambda(members));
    r = numel(members);

    Ealternating = Ealternating + (-1)^(r+1)/L;
end

fprintf('Recursive E[M]      = %.6f hours\n',Erecursive);
fprintf('Alternating E[M]    = %.6f hours\n',Ealternating);
fprintf('Difference          = %.3e\n', ...
    abs(Erecursive-Ealternating));

fprintf('\nExact probability of finishing last:\n');
for i = 1:n
    fprintf('Student %d (rate %.2f): %.6f\n', ...
        i,lambda(i),pLastExact(i));
end

X = -log(rand(N,n))./lambda;
M = max(X,[],2);

[Mmean,~,Mci] = deal(mean(M),[],[]);
Mci = quantile(M,[0.025 0.975]);

fprintf('\nSimulation of %d exams:\n',N);
fprintf('Simulated E[M]       = %.6f hours\n',Mmean);
fprintf('95%% empirical range of M = [%.6f, %.6f]\n', ...
    Mci(1),Mci(2));

[~,firstStudent] = min(X,[],2);
[~,lastStudent] = max(X,[],2);

fprintf('\nFirst-finisher probabilities:\n');
fprintf('Student   Exact       Simulated\n');
for i = 1:n
    fprintf('%3d      %.5f       %.5f\n',i, ...
        lambda(i)/sum(lambda),mean(firstStudent==i));
end

fprintf('\nLast-finisher probabilities:\n');
fprintf('Student   Exact       Simulated\n');
for i = 1:n
    fprintf('%3d      %.5f       %.5f\n',i, ...
        pLastExact(i),mean(lastStudent==i));
end

tMax = max(quantile(M,0.999),max(1./lambda)*5);
t = linspace(0,tMax,1000);

F = ones(size(t));
for i = 1:n
    F = F.*(-expm1(-lambda(i)*t));
end

fM = zeros(size(t));
for i = 1:n
    term = lambda(i)*exp(-lambda(i)*t);
    for j = 1:n
        if j ~= i
            term = term.*(-expm1(-lambda(j)*t));
        end
    end
    fM = fM + term;
end

figure('Color','w');
histogram(M,100,'Normalization','pdf', ...
    'FaceColor',[0.65 0.78 0.92], ...
    'EdgeColor','none');
hold on;
plot(t,fM,'r-','LineWidth',2);
xlabel('Last finishing time M (hours)');
ylabel('Density');
title('Last finishing time: simulation and exact density');
legend('Simulated exams','Exact density','Location','best');
grid on;
hold off;

Xfirst = min(X,[],2);

figure('Color','w');
tiledlayout(2,1);

nexttile;
histogram(Xfirst,80,'Normalization','pdf', ...
    'FaceColor',[0.65 0.78 0.92], ...
    'EdgeColor','none');
hold on;
plot(t, sum(lambda)*exp(-sum(lambda)*t), ...
    'r-','LineWidth',1.8);
xlim([0, min(tMax,0.5)]);
xlabel('First finishing time X');
ylabel('Density');
title('Distribution of the first finishing time');
legend('Simulation','Exact exponential density');
grid on;
hold off;

nexttile;
boxplot(Xfirst,firstStudent);
xlabel('Identity of first finisher');
ylabel('First finishing time X');
title('First finishing time conditional on first finisher');
grid on;

maxK = 12;
Ek = zeros(maxK,1);

for k = 1:maxK
    rates = 1:k;
    total = 2^k;
    Etemp = zeros(total,1);

    for mask = 1:total-1
        members = find(bitget(mask,1:k));
        L = sum(rates(members));

        Etemp(mask+1) = 1/L;
        for i = members
            smallerMask = bitset(mask,i,0);
            Etemp(mask+1) = Etemp(mask+1) + ...
                (rates(i)/L)*Etemp(smallerMask+1);
        end
    end

    Ek(k) = Etemp(total);
end

figure('Color','w');
plot(1:maxK,Ek,'bo-','LineWidth',1.5, ...
    'MarkerFaceColor','b');
hold on;
yline(2*log(2),'r--','2 ln 2','LineWidth',1.5);
xlabel('Number of students k');
ylabel('E[M] (hours)');
title('Expected last finishing time for rates \lambda_i=i');
legend('Exact recursion','Upper bound','Location','best');
grid on;
hold off;

fprintf('\nRates lambda_i = i:\n');
fprintf(' k       E[M]\n');
fprintf('%2d    %.6f\n',[(1:maxK); Ek.']);
fprintf('Upper bound 2 ln 2 = %.6f\n',2*log(2));
