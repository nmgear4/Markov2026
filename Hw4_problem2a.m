p = 0.6;
q = 1 - p;
n = 0:50;

p_return = zeros(size(n));

for k = 1:length(n)
    p_return(k) = nchoosek(2*n(k),n(k)) * (p*q)^n(k);
end

E_N = 1/abs(p-q) - 1;

rho = 1 - abs(p-q);

fprintf('p = %.1f, q = %.1f\n',p,q);
fprintf('E_0[N(0)] = %.2f\n',E_N);
fprintf('rho_00 = %.2f\n',rho);

figure;
plot(2*n,p_return,'o-','LineWidth',1.5);
xlabel('n');
ylabel('p_n(0,0)');
title('Return probabilities for p = 0.6');
grid on;
