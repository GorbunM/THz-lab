omega = linspace(1, 10, 100);
b = omega + sqrt(omega)
c = omega .* sqrt(omega)
for i = 1:100
    k(i,:) = roots([1 -b(i) c(i)])
end

plot(k, omega)