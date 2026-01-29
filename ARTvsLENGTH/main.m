c = 29979245800;
%%
nunu = 3500;
samplMt = dlmread(sprintf('.\\vst_%drcm.txt', nunu));
reflEXP = samplMt(:,2)*1e-2;
reflEXP = reflEXP';
%%
thcEXP   = samplMt(:,1)*1e-6;
thcEXP = thcEXP';
%%
fin = readmatrix('5-thick.txt');
% thcEXP = fin(1, 2:end)*1e-6;
% nstr = length(fin(:,1));
nstr = 269;
Mout = zeros(nstr, 4);
p_opt = [1,0];
for freqidx = 5:nstr
    nunu = fin(freqidx+1, 1);
    reflEXP = fin(freqidx+1, 2:end) * 1e-2;
    %%
    
    err = @(p) sum(log(optima2(p, thcEXP, nunu*c) ./ reflEXP).^2);
    % err = @(p) sum((optima(p, thcEXP) - reflEXP).^2);
    % p0 = [1,0];
    p0 = p_opt;
    % p0 = (p_opt+p0)/2;
    pma = [8, 8];
    pmi = [0, 0];
    [p_opt, fval] = fmincon(err, p0, [],[],[],[], pmi, pma);
    
    
    %%
    
    % ptn=301;
    % thick=linspace(thcEXP(1), thcEXP(end), ptn);
    % mod = optima2(p_opt, thick, nunu*c);
    
    % figure
    % hold on
    % plot(thick, mod)
    % scatter(thcEXP, reflEXP)
    % disp(p_opt)
    
    Mout(freqidx+1, :) = [nunu, p_opt, err(p_opt)];

end
%%
er = Mout(:,4);
epseps = Mout(er < 0.2, :);
figure
% scatter(epseps(:,1), epseps(:,2))
scatter(epseps(:,1), epseps(:,3))

% plot(epseps(:,1), epseps(:,2))
% plot(epseps(:,1), epseps(:,3))

disp('That''s it')