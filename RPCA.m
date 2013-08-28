function [X, Y, funVal] = RPCA(A, z, opts)
% 
% In this function, we approximately solve:
%   min rank(X) + z * ||Y||_0
%   s.t. ||A - X - Y||_F^2 <= \sigma^2
% using alternating method.
% 

% Init.
sigma = opts.sigma;
% this is actually sigma square.

theta1 = opts.theta1;
theta2 = opts.theta2;
% Approximation of 0 norm and rank.

% Init point
switch opts.init
    case 1
        % Starting from 0
        X = A;
        Y = zeros(size(A));
    case 2
        % starting from nuclear + 1 norm optimization
        fprintf('warm start\n');
        tic
        [X, Y] = RobustPCA_sub(A, zeros(size(A)), zeros(size(A)), sigma, 1, 1 / z, zeros(size(A)));
        toc
    case 3
        k = 5;
        [U, G, V] = svd(A, 'econ');
        G(k+1:size(G), :) = 0; 
        X = U*G*V';
        Y = rand(size(A))*sigma;
    case 4
        tol = 1e-5;
        n = size(A, 1);
        dev = sigma/sqrt(n+sqrt(8*n));
        [X, Y]=nsa_v1(A, dev, tol);

end
funVal = zeros(opts.maxIter, 1);

% Start DC programing

for iter = 1:opts.maxIter
    % Check for convergence.
    disp(iter);
    % Fix Y, solve X
    [ X, ~ ] = CalX( theta1, A, Y, sigma );
    % Fix X, solve Y
    [ Y, ~ ] = CalY( theta2, A, X, sigma );
    
    funVal(iter) = RPCA_fun(X, Y);
    disp(funVal(iter));
    if iter >= 2 && abs(funVal(iter) - funVal(iter - 1)) ...
            <= opts.tol * max(1, abs(funVal(iter)))
        funVal(iter + 1:end) = [];
        return;
    end
    
end

    function val = RPCA_fun(X, Y)
        [~, S, ~] = svd(X, 'econ');
        Y = abs(Y(:));
        val = (sum(diag(S)) - sum(max(diag(S) - theta1, 0))) / theta1 ...
            + z * (sum(Y) - sum(max(Y - theta2, 0))) / theta2;
    end

end