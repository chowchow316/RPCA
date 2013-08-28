function [A, X0, Y0, sigma2] = GenSyn_noise(n, cr, cp, noise)
% Without noise
% X0: low rank matrix
% Y0: sparse component
% Z0: noise
% cr: rank coeff for X0
% cp: sparse coeff for Y0
% sigma2: approximation error for constraint

r = round(n * cr);
Range = 100;

U = randn(n, r);
V = randn(n, r);

X0 = U * V';

Y0 = full(sprand(n, n, cp));
Y0(Y0 ~= 0) = - Range + 2 * Range * Y0(Y0 ~= 0);
% Y0 is from [- Range, Range]

Z = noise * randn(n);

sigma2 = (n + sqrt(8 * n)) * noise^2;
sigma2 = sqrt(sigma2);

A = X0 + Y0 + Z;
% 
% disp( norm(Z, 'fro')^2 );
% disp([sigma2, sqrt(sigma2), (n + sqrt(8 * n)) * noise]);