clear, clc;
% Generate A = X + Y;
% 

theta1 = 0.01; %0.01;  
theta2 = 0.01; %0.01

addpath(genpath('./exact_alm_rpca'));
addpath(genpath('./NSA_v2'));

%a = [0.05, 0.1];
n = 100; % 100, 200, 500.

cr = 0.1; % 0.05
cp = 0.05;

noise = 0.001; %1e-3

[A, X0, Y0, sigma2] = GenSyn_noise(n, cr, cp, noise);

z = 1 / sqrt(n);
tol = 1e-5;

% our method

opts = [];
opts.theta1 = theta1;
opts.theta2 = theta2;
opts.sigma = sigma2;
opts.sub = 'ADMM';
opts.init = 4;  %4 for closed form
opts.tol = tol;
opts.maxIter = 500;

% ASALM
tic
[X1, Y1] = exact_alm_rpca(A, z, tol);
toc


% NSA
tic
%[X2, Y2] = NSA(A, sigma2, z);
dev = sigma2/sqrt(n+sqrt(8*n));
[X2, Y2]=nsa_v1(A, dev, tol);
toc


% closed form
tic;
[X, Y, funVal] = RPCA(A, z, opts);
toc;

opts.init = 2;  %4 for closed form
% DC programming
% tic
% [X3, Y3, funVal1] = RobustPCA(A, z, opts);
% toc
% disp(norm(X3 + Y3 - A, 'fro')/norm(A, 'fro'));
% disp(norm(X3 - X0, 'fro') / norm(X0, 'fro'));
% disp(norm(Y3 - Y0, 'fro') / norm(Y0, 'fro'));
% disp(nnz((Y3~= 0) == (Y0 ~= 0)) / n^2);
% disp([rank(X3),  nnz(Y3)]);



% disp results
% ALM
disp('ALM');
disp(norm(X1 + Y1 - A, 'fro')/norm(A, 'fro'));
disp(norm(X1 - X0, 'fro') / norm(X0, 'fro'));
disp(norm(Y1 - Y0, 'fro') / norm(Y0, 'fro'));
disp(nnz((Y1 ~= 0) == (Y0 ~= 0)) / n^2);
disp([rank(X1), nnz(Y1)]);

%  NSA
disp('NSA');
disp(norm(X2 + Y2 - A, 'fro')/norm(A, 'fro'));
disp(norm(X2 - X0, 'fro') / norm(X0, 'fro'));
disp(norm(Y2 - Y0, 'fro') / norm(Y0, 'fro'));
disp(nnz((Y2 ~= 0) == (Y0 ~= 0)) / n^2);
disp([rank(X2), nnz(Y2)]);


% OUR
disp('Our');
disp(norm(X + Y - A, 'fro')/norm(A, 'fro'));
disp(norm(X - X0, 'fro') / norm(X0, 'fro'));
disp(norm(Y - Y0, 'fro') / norm(Y0, 'fro'));
disp(nnz((Y ~= 0) == (Y0 ~= 0)) / n^2);
disp([rank(X),  nnz(Y)]);

% [~, S, ~] = svd(X);
% diag_S = diag(S);
% 
% Ya = abs(Y(:));
% disp((sum(diag_S) - sum(max(diag_S - theta1, 0))) / theta1 ...
%     + z * (sum(Ya) - sum(max(Ya - theta2, 0))) / theta2);
% disp(rank(X) + z * nnz(Y));
% 
% [~, S, ~] = svd(X1);
% diag_S = diag(S);
% 
% Y1a = abs(Y(:));
% disp((sum(diag_S) - sum(max(diag_S - theta1, 0))) / theta1 ...
%     + z * (sum(Y1a) - sum(max(Y1a - theta2, 0))) / theta2);
% disp(rank(X1) + z * nnz(Y1));
% 
% [~, S, ~] = svd(X2);
% diag_S = diag(S);
% 
% Y2a = abs(Y2(:));
% disp((sum(diag_S) - sum(max(diag_S - theta1, 0))) / theta1 ...
%     + z * (sum(Y2a) - sum(max(Y2a - theta2, 0))) / theta2);
% disp(rank(X2) + z * nnz(Y2));
 %   end
%end
%save results/11synthetic100.mat A sigma2 X0 Y0 X Y X1 Y1 X2 Y2;