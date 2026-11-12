Q = [-2  2  0  0;
      1 -3  2  0;
      0  2 -4  2;
      0  0  3 -3];

lambda = sort(eig(Q), 'descend');

disp('Eigenvalues:');
disp(lambda);

disp('Sum of eigenvalues:');
disp(sum(lambda));

disp('Trace of Q:');
disp(trace(Q));

lambda2 = lambda(2);

pi = [3; 6; 6; 4]/19;

disp('Stationary distribution:');
disp(pi);

x0 = [1; 0; 0; 0];

t = linspace(0, 5, 1001);

X = zeros(4, length(t));

TV = zeros(1, length(t));

for j = 1:length(t)

    X(:,j) = expm(Q' * t(j)) * x0;

    TV(j) = 0.5 * sum(abs(X(:,j) - pi));

end

figure;

semilogy(t, X(1,:), 'LineWidth', 2);
hold on;
semilogy(t, X(2,:), 'LineWidth', 2);
semilogy(t, X(3,:), 'LineWidth', 2);
semilogy(t, X(4,:), 'LineWidth', 2);

xlabel('Time (hours)');
ylabel('Probability (log scale)');

title('State Probabilities Over Time');

legend('State 0', 'State 1', ...
       'State 2', 'State 3', ...
       'Location', 'best');

grid on;

figure;

semilogy(t, TV, 'b', 'LineWidth', 2);
hold on;

reference = TV(1) * exp(lambda2 * t);

semilogy(t, reference, 'r--', 'LineWidth', 2);

xlabel('Time (hours)');
ylabel('TV distance (log scale)');

title('Convergence to Stationary Distribution');

legend('TV(t)', ...
       'Reference slope \lambda_2', ...
       'Location', 'best');

grid on;

TVfun = @(s) 0.5 * sum(abs(expm(Q' * s) * x0 - pi));

t01 = fzero(@(s) TVfun(s) - 0.01, [0 5]);

fprintf('Time when TV = 0.01: %.6f hours\n', t01);

t_approx = log(100)/abs(lambda2);

fprintf('Approximation ln(100)/|lambda2|: %.6f hours\n', ...
        t_approx);

fprintf('Difference: %.6f hours\n', t_approx - t01);
