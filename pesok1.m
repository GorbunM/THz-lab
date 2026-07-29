

c = 29979245800; % speed of light in cm/s 
omega = 1e12;

sigma = 2e14;
e1 = 1;
e2 = 2;

s = sigma/c;

Coef(1) = 256*pi^4*s^4;

Coef(2) = -32*(16*pi^4*e1*s^4 + 16*pi^4*e2*s^4 - pi^2*e1^2*s^2 - pi^2*e2^2*s^2);

Coef(3) = 256*pi^4*e1^2*s^4 + 1024*pi^4*e1*e2*s^4 + 256*pi^4*e2^2*s^4 - 32*pi^2*e1^3*s^2 ...
    - 64*pi^2*e1^2*e2*s^2 - 64*pi^2*e1*e2^2*s^2 - 32*pi^2*e2^3*s^2 + e1^4 - 2*e1^2*e2^2 + e2^4;

Coef(4) = 2*(256*pi^4*e1^2*e2*s^4 + 256*pi^4*e1*e2^2*s^4 - 32*pi^2*e1^3*e2*s^2 ...
    - 32*pi^2*e1^2*e2^2*s^2 - 32*pi^2*e1*e2^3*s^2 + e1^4*e2 - e1^3*e2^2 - e1^2*e2^3 + e1*e2^4);

Coef(5) = 256*pi^4*e1^2*e2^2*s^4 - 32*pi^2*e1^3*e2^2*s^2 - 32*pi^2*e1^2*e2^3*s^2 + e1^4*e2^2 ...
    - 2*e1^3*e2^3 + e1^2*e2^4;

for i = 1
    a(i,:) = roots(Coef)
end

k_p = sqrt(a)*omega/c