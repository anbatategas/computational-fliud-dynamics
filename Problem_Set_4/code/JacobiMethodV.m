function [x_final, iter] = JacobiMethodV(A, b, x_initial, epsilon)
% Function JacobiMethodV.m
%---------------------------------------------------------------
% Jacobi's iterative method for linear system A*x = b in vector form.
% The vector form is used for reducing computation cost and increasing
% computation speed.
% Input : - square matrix A of coefficients   (n x n).
%         - const. vector b                   (n x 1).
%         - inital guess x_initial            (n x 1).
%         - termination parameter epsilon     (1 x 1).
%           (use of "weighted" infinity norm)
% Output: - vector x_final = sol. of (Ax=b)   (n x 1).
%         - integer iter = no. of iterations  (1 x 1).
%---------------------------------------------------------------

% --- System dimension & solution definition ---
n = length(b);

% --- Input validity ---
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
diag_A = diag(A);
sum_rows = sum(abs(A), 2) - abs(diag_A);
if any(sum_rows > abs(diag_A))
    warning('Matrix A is not strictly diagonally dominant. Method may diverge.');
end
if any(diag_A == 0)
    error('Zero diagonal element detected.');
end
% Ensure A has nonzero diagonal elements:
if(any(diag(A) == 0))
    error('Matrix A has a zero as diagonal element. Method diverges');
end

% --- Jacobi's method implementation ---
% --- 0. Matrix separation into D and R ---
D = diag_A;              % Vector of diagonal elements
R = A - diag(D);         % Matrix A with zeros on the diagonal
D_inv = 1 ./ D;          % Inverse of diagonal (element-wise division)
% --- 1. Temporary and Input variables initialization ---
x_old = x_initial;
x_final = x_old;
iter = 0;
terminationRatio = 1 + epsilon;
% --- 2. Main iteration loop ---
while(terminationRatio >= epsilon)
    x_final = D_inv .* (b - R * x_old);
    iter = iter + 1;
    terminationRatio = norm(x_final - x_old, inf) / norm(x_final, inf);
    x_old = x_final;
end
end