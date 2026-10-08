clear
clc

% --- CFD Exercise 1.3 ---
% Solution script (using NRMethod.m)
%
% --- Definition of function handle ---
F = @(x) [x(1)^2+x(2)-37;
          x(1)-x(2)^2-5;
          x(1)+x(2)+x(3)-3]
%
% --- Definition of x_initial ---
x_initial = [0 0 0]
%
% --- Definition of accuracy epsilon ---
epsilon = 1e-5
% --- Apply Newton-Raphson method and get solution ---
[x,iter] = NRMethod(F, x_initial, epsilon)
% End of script.