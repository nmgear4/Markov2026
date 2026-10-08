a = 1;
b = 1;
m = 4;
N = 60;

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

q0 = [1 0 0 0 0];

Q = zeros(N+1,m+1);
Q(1,:) = q0;

for n = 1:N
    Q(n+1,:) = Q(n,:)*P;
end

fprintf('\nq_50 = \n');
disp(Q(51,:));

fprintf('q_51 = \n');
disp(Q(52,:));

runningAvg = zeros(N+1,1);

for n = 0:N
    runningAvg(n+1) = mean(Q(1:n+1,3));
end

n = 0:N;

figure;

subplot(3,1,1);
plot(n,Q(:,3),'o-','LineWidth',1.2);
xlabel('n');
ylabel('q_n(2)');
title('q_n(2)');
grid on;

subplot(3,1,2);
plot(n,Q(:,5),'o-','LineWidth',1.2);
xlabel('n');
ylabel('q_n(4)');
title('q_n(4)');
grid on;

subplot(3,1,3);
plot(n,runningAvg,'LineWidth',1.5);
xlabel('n');
ylabel('Running average');
title('Running average of q_n(2)');
grid on;
