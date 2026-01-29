thcEXP = thclist*1e-2;
P1 = zeros(npt, 2);
Err = zeros(npt, 1);
% Ropt = zeros(npt, )
for nunu = 1:npt
    reflEXP = Rout(:,nunu)';
    nu0 = nu(nunu);
    
    err = @(p) sum(log(optima2(p, thcEXP, nu0) ./ reflEXP).^2);
    % err = @(p) sum((optima(p, thcEXP) - reflEXP).^2);
    p0 = [3.6,0];
    % p0 = p_opt;
    % p0 = p1;
    pma = [10, 10]*100;
    pmi = pma*(-1);
    % pmi = [3,0];
    [p_opt, fval] = fmincon(err, p0, [], [], [], [], pmi, pma);
    Err(nunu) = fval;
    P1(nunu, :) = p_opt;
end

disp('ready')
%%

Ropt = zeros(nthc, npt);
P1(:, 2) = abs(P1(:, 2));
for nunu = 1:npt
    reflEXP = Rout(:,nunu)';
    nu0 = nu(nunu);
    Ropt(:, nunu) = optima2(P1(nunu,:), thcEXP, nu0);
end

disp('ready')
%%

figure
surf(nu, thclist, Rout, 'EdgeColor','none')
set(gca,'YScale','log');
xlabel('\nu');
ylabel('d');
zlabel('R');
figure

surf(nu, thclist, Ropt, 'EdgeColor','none')
set(gca,'YScale','log');
xlabel('\nu');
ylabel('d');
zlabel('R');

%%
figure
hold on
plot(real(eps_layers(2,:)))
plot(P1(:,1))

figure
hold on
plot(imag(eps_layers(2,:)))
plot(P1(:,2))