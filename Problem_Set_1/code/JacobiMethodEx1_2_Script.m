clear
clc

% --- CFD Exercise 1.2 ---
% Solution script (using JacobiMethod.m)
%
% --- Definition of matrix A ---
A = [10 5 0 0;
    5 10 -4 0;
    0 -4 8 -1;
    0 0 -1 5]
%
% --- Definition of constant vector b ---
b = [6 25 -11 -11]
% Note: It is ok to define b as a row. The function
% in "JacobiMethod.m" will take care of it.
%
% --- Definition of x_initial ---
x_initial = [0 0 0 0]
%
% --- Definition of accuracy epsilon ---
epsilon = 1e-5
% --- Apply Jacobi method and get solution ---
[x,iter] = JacobiMethod(A,b,x_initial,epsilon)
%
% --- Check validity of solution ---
x1 = A\b'
if(abs(x-x1)<epsilon)
    disp(['Solution via Jacobi method and via A\b match ' ...
        'within epsilon chosen!'])
else
    disp(['Solution via Jacobi method and via A\b differ ' ...
        'at most by epsilon chosen!'])
end
% End of script.