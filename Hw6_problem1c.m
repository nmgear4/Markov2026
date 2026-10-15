fprintf('%3s %6s %16s %14s\n', ...
    'm', 'M', 'Wielandt bound', '|lambda_2|');

for m = 3:12

    P = zeros(m,m);

    for i = 1:m-1
        P(i,i+1) = 1;
    end

    P(m,1) = 1/2;
    P(m,2) = 1/2;

    A = double(P > 0);
    B = A;
    M = 1;

    while ~all(B(:))
        B = double((B * A) > 0);
        M = M + 1;
    end

    lambda = eig(P);

    lambda = lambda(abs(lambda - 1) > 1e-7);

    lambda2 = max(abs(lambda));

    bound = (m-1)^2 + 1;

    fprintf('%3d %6d %16d %14.6f\n', ...
        m, M, bound, lambda2);
end
