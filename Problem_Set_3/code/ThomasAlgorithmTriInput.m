function [x] = ThomasAlgorithmTriInput(a,b,c,d)
% Function ThomasAlgorithmTriInput.m
%---------------------------------------------------------------
% Solution of tridiagonal linear system A*x = d, A=diag(a,b,c)
% Input : - diagonal vectors a,b,c of matrix A          (n x 1).
%         - const. vector d                             (n x 1).
% Output: - vector x = solution of (Ax=d)               (n x 1).
%---------------------------------------------------------------

% --- System dimension & solution definition ---
n = length(d);
x = zeros(n,1);

% --- Input validity ---
% Ensure d is in column form :
d = d(:);
% Ensure a,b,c and d are properly sized :
if(length(a) ~= n || a(1) ~= 0)
    error('Lower diagonal a is not (n x 1) or a(1) isnot 0.');
elseif(length(b) ~= n)
    error('Diagonal b is not (n x 1).');
elseif(length(c) ~= n || c(n) ~= 0)
    error('Upper diagonal c is not (n x 1) or c(n) isnot 0.');
end

% --- Thomas' method implementation ---
% --- 1. Forward Elimination Stage ---
c(1) = c(1)/b(1);
d(1) = d(1)/b(1);
% Loop over vector c:
for i = 2 : n-1
    denominator = b(i)-a(i)*c(i-1);
    if(denominator==0)
       error('Matrix A is singular or has bad condition.'); 
    end
    c(i) = c(i)/denominator;
end
% Loop over vector d:
for i = 2 : n
    denominator = b(i)-a(i)*c(i-1);
    if(denominator==0)
       error('Matrix A is singular or has bad condition.'); 
    end
    d(i) = (d(i)-a(i)*d(i-1))/denominator;
end
x(n) = d(n);
% --- 2. Backward Substitution Stage ---
for i = n-1 : -1 : 1
    x(i) = d(i) - c(i)*x(i+1);
end
end