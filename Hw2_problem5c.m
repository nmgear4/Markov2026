N = 1e4;
R = 1e3;

C = ones(R,1);

for n =3:(N-1)
        p = C/n;
        U = rand(R,1);
        C = C + (U < p);
end

z = C/N;

emp_m = mean(z);
emp_sd = std(z);
emp_ratio = emp_sd/emp_m;

fprintf('Empirical mean = %.5f\n', emp_m);
fprintf('Empirical SD/mean = %.5f\n', emp_ratio);

minc = min(C);
maxc = max(C);

fprintf('Smallest core = %d\n', minc);
fprintf('Largest core = %d\n', maxc);

binWidth = 0.02;

histogram (z, 'BinWidth', binWidth, 'Normalization', 'pdf');
hold on

zz = linspace(0,1,500);
plot(zz, 2*(1-zz), 'r', 'LineWidth', 2);

xlabel('z=C/N');
ylabel('Density');
title('Distribution of C/N at N = 10^4');
legend('Empirical PDF', 'h(z) = 2(1-z)');
grid on
