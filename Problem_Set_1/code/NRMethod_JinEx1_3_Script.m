clear
clc

% --- CFD Exercise 1.3 ---
% Solution script (using NRMethod.m , NRMethod_Jin.m)
%
% --- Definition of function handle F ---
F = @(x) [x(1)^2+x(2)-37;
          x(1)-x(2)^2-5;
          x(1)+x(2)+x(3)-3]
%
% --- Definition of function handle J ---
J = @(x) [2*x(1),1,0; 1,-2*x(2),0; 1,1,1]

% --- Definition of x_initial ---
x_initial = [0 0 0]
%
% --- Definition of accuracy epsilon ---
epsilon = 1e-5
% --- Apply Newton-Raphson method and get solution ---
[xNR,iterNR] = NRMethod(F, x_initial, epsilon)
% --- Apply Newton-Raphson (with J) method and get solution ---
[xNRJ,iterNRJ] = NRMethod_Jin(F, J, x_initial, epsilon)
% End of script.