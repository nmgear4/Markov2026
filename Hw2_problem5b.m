U = rand(1e5,1);
Z = 1 - sqrt(1-U);

fprintf('Sample mean = %.5f\n', mean(Z));
