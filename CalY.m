function [ Y, objY ] = CalY( theta2, A, X, sigma2 )
% Fix X, calculate Y
% use closed form
[m, n] = size(A);
Z = A - X;
% vectorize
Z = Z(:);
flag = sign(Z);
% abs
Z = abs(Z);
% sort
[Z, ind] = sort(Z); 
% ref = (1: 1: n*n)';
Y = Z;
i = 1;

while sigma2 > 0
    if i > m*n
        break;
    end    
    if Z(i)^2 <= sigma2
        Y(i) = 0;
        sigma2 = sigma2 - Z(i)^2; 
        i = i + 1;
    else
        Y(i) = Z(i) - sqrt(sigma2);
        sigma2 = 0;
    end    

end

objY = sum(min(Y, theta2));
Y(ind) = Y;
% [~, index] = sort(ind);
% id = ref(index);
% Y = Y(id);
Y = flag.*Y;
Y = vec2mat(Y, m)';
    
end

