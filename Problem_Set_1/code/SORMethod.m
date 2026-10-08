function [x_final, iter] = SORMethod(A, b, x_initial, epsilon, omega)
% Function SORMethod.m
%---------------------------------------------------------------
% SOR iterative method for linear system A*x = b.
% Input : - square matrix A of coefficients   (n x n).
%         - const. vector b                   (n x 1).
%         - inital guess x_initial            (n x 1).
%         - termination parameter epsilon     (1 x 1).
%           (use of "weighted" infinity norm)
%         - relaxation parameter omega        (1 x 1).
% Output: - vector x_final = sol. of (Ax=b)   (n x 1).
%         - integer iter = no. of iterations  (1 x 1).
%---------------------------------------------------------------

% --- System dimension & solution definition ---
n = length(b);

% --- Input validity ---
% Ensure omega is between 0 and 2:
if(omega <= 0 || omega >=2)
    error('Relaxation parameter omega must be in interval (0,2).');
end
% Ensure epsilon is positive:
if(epsilon <= 0)
    error('Termination parameter epsilon must be positive.');
end
% Ensure b and x_initial are in column form :
b = b(:);
x_initial = x_initial(:);
% Ensure A and b are properly sized :
[m1,m2]=size(A);
if(m1 ~= m2)
    error('Matrix A is not square (n x n).');
elseif(m1 ~= n)
    error('Matrix A and vector b are incompatible.');
elseif(m1 ~= length(x_initial))
    error('Initial guess x_initial is incompatible with matrices A,b.');
end
% Ensure A is diagonally dominant:
Asumrow = zeros(n,1);
counter = 0;
for i = 1:n
    Asumrow(i) = sum(abs(A(i,:))) - abs(A(i,i));
    if(Asumrow(i) <= abs(A(i,i)))
        counter = counter + 1;
    end
end
if(counter ~= n)
    error('Matrix A is not diagonally dominant. Method may diverge.');
end
% Ensure A has nonzero diagonal elements:
if(any(diag(A) == 0))
    error('Matrix A has a zero as diagonal element. Method diverges');
end

% --- SOR method implementation ---
% --- 1. Temporary and Input variables initialization ---
x_old = x_initial;
x_final = x_old;
iter = 0;
terminationRatio = 1 + epsilon;
% --- 2. Main iteration loop ---
while(terminationRatio >= epsilon)
    for i = 1:n
        s = 0;
        for j = 1:n
            if(j < i)
                s = s + A(i,j)*x_final(j);
            elseif(j > i)
                s = s + A(i,j)*x_old(j);
            end
        end
        x_final(i) = omega*(b(i) - s)/A(i,i) + (1-omega)*x_old(i);
    end
    iter = iter + 1;
    terminationRatio = norm(x_final-x_old, "inf")/norm(x_final,"inf");
    x_old = x_final;
end
end