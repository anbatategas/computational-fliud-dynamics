clear
clc

% --- CFD Exercise 1.2 ---
% Solution script (using JacobiMethod.m, GaussSeidelMethod.m,
%                                                  SORMethod.m)
%
% --- Definition of matrix A ---
A = [10 5 0 0;
    5 10 -4 0;
    0 -4 8 -1;
    0 0 -1 5]
%
% --- Definition of constant vector b ---
b = [6 25 -11 -11]
% Note: It is ok to define b as a row. The functions
% used will take care of it.
%
% --- Definition of x_initial ---
x_initial = [0 0 0 0]
%
% --- Definition of accuracy epsilon ---
epsilon = 1e-5
% --- Apply Jacobi method and get solution ---
[xJ,iterJ] = JacobiMethod(A,b,x_initial,epsilon)
%
% --- Apply Gauss-Seidel method and get solution ---
[xGS,iterGS] = GaussSeidelMethod(A,b,x_initial,epsilon)
%
% --- Apply SOR method (omega=1.1) and get solution ---
[xSOR,iterSOR] = SORMethod(A,b,x_initial,epsilon,1.1)
%
% --- Check validity of solutions ---
x1 = A\b'
if(abs(xJ-x1)<epsilon)
    disp(['Solution via Jacobi method and via A\b match ' ...
        'within epsilon chosen!'])
else
    disp(['Solution via Jacobi method and via A\b differ ' ...
        'at most by epsilon chosen!'])
end
if(abs(xGS-x1)<epsilon)
    disp(['Solution via Gauss-Seidel method and via A\b match ' ...
        'within epsilon chosen!'])
else
    disp(['Solution via Gauss-Seidel method and via A\b differ ' ...
        'at most by epsilon chosen!'])
end
if(abs(xSOR-x1)<epsilon)
    disp(['Solution via SOR method and via A\b match ' ...
        'within epsilon chosen!'])
else
    disp(['Solution via SOR method and via A\b differ ' ...
        'at most by epsilon chosen!'])
end
% End of script.