R0vals = linspace(0.5, 3, 500);
Nvals = [100, 1000];

figure('Color','w');
hold on;

colors = lines(3);

for a = 1:length(Nvals)
    N = Nvals(a);
    K = N/10;
    h1 = zeros(size(R0vals));

    for i = 1:length(R0vals)
        R0 = R0vals(i);

        rho = ones(K,1);

        for j = 1:K-1
            rho(j+1) = rho(j) * ...
                N/(R0*(N-j));
        end

        h1(i) = 1/sum(rho);
    end

    plot(R0vals,h1,'LineWidth',2,...
        'Color',colors(a,:),...
        'DisplayName',sprintf('Exact, N=%d',N));
end

hLimit = max(0,1 - 1./R0vals);

plot(R0vals,hLimit,'k--','LineWidth',2,...
    'DisplayName','Limit: max(0,1-1/R_0)');

xlabel('R_0');
ylabel('Probability of reaching K before extinction');
title('Outbreak probability from one initial infection');
legend('Location','northwest');
grid on;
xlim([0.5 3]);
ylim([0 0.75]);

exportgraphics(gcf,'outbreak_probability.png',...
    'Resolution',200);
