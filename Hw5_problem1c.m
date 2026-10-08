a = 0.3;
b = 0.1;
m = 4;
N = 200;

P = zeros(m+1);

for k = 0:m
    i = k + 1;

    if k < m
        P(i,i+1) = a*(m-k)/m;
    end

    if k > 0
        P(i,i-1) = b*k/m;
    end

    P(i,i) = 1 - sum(P(i,:));
end

disp('Transition matrix P =');
disp(P);

A = P' - eye(m+1);
A(end,:) = ones(1,m+1);

rhs = zeros(m+1,1);
rhs(end) = 1;

pi = A\rhs;

fprintf('\nStationary distribution from linear solve:\n');
disp(pi');

theta = a/(a+b);

pi_binomial = zeros(1,m+1);

for k = 0:m
    pi_binomial(k+1) = nchoosek(m,k) * ...
        theta^k * (1-theta)^(m-k);
end

fprintf('Stationary distribution from binomial formula:\n');
disp(pi_binomial);

fprintf('Maximum difference between the two:\n');
disp(max(abs(pi' - pi_binomial)));

q0 = [1 0 0 0 0];

Q = zeros(N+1,m+1);
Q(1,:) = q0;

for n = 1:N
    Q(n+1,:) = Q(n,:)*P;
end

tol = 1e-6;
n_min = NaN;

for n = 0:N
    err = max(abs(Q(n+1,:) - pi'));
    
    if err < tol
        n_min = n;
        break;
    end
end

fprintf('\nSmallest n satisfying the tolerance:\n');
fprintf('n = %d\n',n_min);


eigenvalues = eig(P);

[~,idx] = sort(abs(eigenvalues),'descend');
eigenvalues_sorted = eigenvalues(idx);

SLEM = abs(eigenvalues_sorted(2));

fprintf('\nEigenvalues of P:\n');
disp(eigenvalues);

fprintf('Second-largest eigenvalue modulus (SLEM): %.6f\n',SLEM);


n = 0:N;
