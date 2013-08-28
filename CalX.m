function [ X, objX ] = CalX( theta1, A, Y, sigma2 )
% Fix Y, solve X
% use closed form
% [m, n] = size(A);
Z = A - Y;
[U, S, V] = svd(Z, 'econ'); % econ SVD
% vectorize
Z = diag(S);
% sort
[Z, ind] = sort(Z);
% ref = (1: 1: length(Z))';
X = Z;
i = 1;
while sigma2 > 0
    if i > length(Z)
        break;
    end
    if Z(i)^2 <= sigma2
        X(i) = 0;
        sigma2 = sigma2 - Z(i)^2; 
        i = i + 1;
    else
        X(i) = Z(i) - sqrt(sigma2);
        sigma2 = 0;
    end
    
end
objX = sum(min(X, theta1));
% [~, index] = sort(ind);
% id = ref(index);
% X = X(id);
X(ind) = X;
X = diag(X);
X = U*X*V';

end

