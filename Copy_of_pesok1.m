A = [];
for k = 1:200
    [a, b, c, Abs, d, e] = parametricTMM(X(k));
    A = [A; Abs]; 
end

figure
plot(A)