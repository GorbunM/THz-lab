err = @(p) sum(log(optima2(p, thcEXP, nu) ./ reflEXP).^2);

p0 = [1,0];
pma = [10, 10];
pmi = [0, 0];
[p_opt, fval] = fmincon(err, p0, [],[],[],[], pmi, pma);

disp(p_opt)