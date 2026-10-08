clear;
clc;
% --- CFD Exercise 2.1 ---
%      --- Part 3 ---
% Solution script (using RK2Method.m)
%
% --------------------------------
% --- DEFINITION OF PARAMETERS ---
% --------------------------------
% --- Definition of function handle ---
f = @(t,y) 1 + (t-y).^2;
% --- Definition of t_initial, t_final ---
t_initial = 2; t_final = 3;
% --- Definition of y0 ---
y0 = 1;
% --- Definition of partition step h ---
h = 0.1;
% --- Definition of accuracy epsilon ---
epsilon = 1e-5;
% --- Definition of number N of intervals ---
N = floor((t_final-t_initial)/h);
% --- Definition of node number vector ---
i = (1:(N+1))';
% --- Definition of time values vector ---
% (this vector will be overwritten below, but it is the same
%                         for given h, t_initial and t_final)
t = (t_initial : h : t_final)';
% --- Definition of analytical solution function handle ---
y_analytical = @(t) t + (1-t).^-1;
% --- Definition of analytical solution at partition nodes ---
y_anal = y_analytical(t);
%
% --------------------------
% --- SOLUTION OF I.V.P. ---
% --------------------------
% --- Runge-Kutta 2nd order method solution ---
[t,y_RK2] = RK2Method(f,t_initial,t_final,y0,h);
% --- Error Calculation ---
error_RK2 = abs(y_RK2 - y_anal);
% --- Definition of table for data output ---
T = table(i, t, y_anal, y_RK2, error_RK2);
disp(T);
% --- Plot solution and analytical values vs. t values ---
figure;
plot(t, y_anal, 'b-', 'LineWidth', 0.5); hold on;
plot(t, y_RK2, 'rx', 'LineWidth', 0.5);
grid on;
xlabel('Time (t)');
ylabel('Solution y(t)');
yticks(1: 0.2 : 2.6);
title('Runge-Kutta 2nd Order method (h = 0.1)');
legend('Analytical Solution', 'R-K 2nd Order Method', 'Location', 'Best');
% End of script.