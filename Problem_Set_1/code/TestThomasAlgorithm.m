function [x] = TestThomasAlgorithm(A,d) 
% Function ThomasAlgorithm.m 
%--------------------------------------------------------------- 
% Solution of tridiagonal linear system A*x = d. 
% Input : - square tridiagonal matrix A of coefficients (n x n). 
%         - const. vector d                             
% Output: - vector x = solution of (Ax=d)               
%--------------------------------------------------------------- 
% --- System dimension & solution definition --- 
n = length(d); 
x = zeros(n,1); 
% --- Input validity --- 
% Ensure d is in column form : 
d = d(:); 
% Ensure A and d are properly sized : 
[m1,m2]=size(A); 
if(m1 ~= m2) 
error('Matrix A is not square (n x n).'); 
elseif(m1 ~= n) 
error('Matrix A and vector d are incompatible.'); 
end 
% Ensure A is a tridiagonal matrix: 
if(~isbanded(A,1,1)) 
error('Matrix A is not tridiagonal.'); 
end 
% --- Extract from A only the non-zero diagonals --- 
% Form vector a (n x 1) of lower-diagonal elements: 
a = zeros(n,1); 
a(2:n) = diag(A,-1); 
% Form vector b (n x 1) of diagonal elements: 
b = zeros(n,1); 
b = diag(A,0); 
% Form vector c (n x 1) of upper-diagonal elements: 
c = zeros(n,1); 
c(1:n-1) = diag(A,1); 
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