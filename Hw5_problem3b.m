P = [ ...
    0    1/2  1/2  0    0    0    0    0;
    0    0    1/2  1/2  0    0    0    0;
    1/2  0    0    0    1/2  0    0    0;
    1/3  0    1/3  0    1/3  0    0    0;
    0    1/3  0    0    0    1/3  1/3  0;
    0    0    0    0    0    1    0    0;
    0    0    0    0    0    0    0    1;
    0    0    0    0    0    0    1    0];

d = 0.85;
n = 8;

G = d*P + (1-d)/n * ones(n);

fprintf('Google matrix G =\n');
disp(G);

q = ones(1,n)/n;

tol = 1e-10;
iteration = 0;

while true

    qnew = q*G;
    iteration = iteration + 1;

    error = norm(qnew-q,1);

    if error < tol
        break;
    end

    q = qnew;
end

pi = qnew;

fprintf('\nNumber of power iterations = %d\n',iteration);

fprintf('\nPageRank vector from power iteration:\n');
disp(pi);

fprintf('Sum of PageRank = %.12f\n',sum(pi));

fprintf('\nCheck pi*G - pi:\n');
disp(pi*G - pi);

A = G' - eye(n);

A(n,:) = ones(1,n);

rhs = zeros(n,1);
rhs(n) = 1;

pi_linear = (A\rhs)';

fprintf('\nPageRank from linear solve:\n');
disp(pi_linear);

fprintf('\nDifference between power iteration and linear solve:\n');
disp(pi - pi_linear);

[sortedValues,ranking] = sort(pi,'descend');

fprintf('\nRanking from highest to lowest:\n');

for k = 1:n
    fprintf('%d. Page %d   PageRank = %.8f\n', ...
        k, ranking(k), sortedValues(k));
end
