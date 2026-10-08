clear; clc;
%      --- CFD Exercise 2.1 ---
% --- NUMERICAL METHODS COMPARISON ---
%
% (using EulerExplicitMethod.m,
%        EulerImplicitMethod.m,
%        RK2Method.m, RK3Method.m, RK4Method.m)
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
% ---------------------------
% ---  SOLUTION OF I.V.P. ---
% --- (USING ALL METHODS) ---
% ---------------------------
% --- 1. Euler's Explicit method solution ---
[~,y_EE] = EulerExplicitMethod(f,t_initial,t_final,y0,h);
% --- Error Calculation ---
error_EE = abs(y_EE - y_anal);
% --- 2. Euler's Implicit method solution ---
[~,y_EI] = EulerImplicitMethod(f,t_initial,t_final,y0,h,epsilon);
% --- Error Calculation ---
error_EI = abs(y_EI - y_anal);
% --- 3. Runge-Kutta 2nd order method solution ---
[~,y_RK2] = RK2Method(f,t_initial,t_final,y0,h);
% --- Error Calculation ---
error_RK2 = abs(y_RK2 - y_anal);
% --- 4. Runge-Kutta 3rd order method solution ---
[~,y_RK3] = RK3Method(f,t_initial,t_final,y0,h);
% --- Error Calculation ---
error_RK3 = abs(y_RK3 - y_anal);
% --- 5. Runge-Kutta 4th order method solution ---
[~,y_RK4] = RK4Method(f,t_initial,t_final,y0,h);
% --- Error Calculation ---
error_RK4 = abs(y_RK4 - y_anal);
% -------------------------------------------
% --- Definition of table for data output ---
l = repmat('|',N+1,1);
T1 = table(i,l, t,l, y_anal,l, y_EE,error_EE,l, y_EI,error_EI,l);
T2 = table(l, y_RK2,error_RK2,l, y_RK3,error_RK3,l, y_RK4,error_RK4);
T1 = renamevars(T1,"l"," "); T1 = renamevars(T1,"l_1","  ");
T1 = renamevars(T1,"l_2","   "); T1 = renamevars(T1,"l_3","    ");
T1 = renamevars(T1,"l_4"," ... "); T2 = renamevars(T2,"l"," ... ");
T2 = renamevars(T2,"l_1","    "); T2 = renamevars(T2,"l_2","     ");
disp(T1); disp(T2);
% -----------------------------------------------
% --- ERROR COMPARISON PLOTS ---
% -----------------------------------------------
% --- Plot solution error values vs. t values ---
figure;
semilogy(t, error_EE,  'r+', 'LineWidth', 0.5); hold on;
semilogy(t, error_EI,  'go', 'LineWidth', 0.5); hold on;
semilogy(t, error_RK2, 'bx', 'LineWidth', 0.5); hold on;
semilogy(t, error_RK3, 'm^', 'LineWidth', 0.5); hold on;
semilogy(t, error_RK4, 'kv', 'LineWidth', 0.5);
grid on;
xlabel('Time (t)');
ylabel('Absolute Error');
title('Numerical Method Comparison (h = 0.1)');
legend('Euler''s Explicit Method', ...
    'Euler''s Implicit Method','R-K 2nd Order Method', ...
    'R-K 3rd Order Method (Rev)', 'R-K 4th Order Method', ...
    'Position', [0.55,0.25,0.35,0.18]);
% ------------------------------------------------
% --- Plot max abs error value vs. step size h ---
h = [0.1 0.05 0.025 0.01 0.005 0.0025 0.001];
max_err_EE  = zeros(length(h),1);
max_err_EI  = zeros(length(h),1);
max_err_RK2 = zeros(length(h),1);
max_err_RK3 = zeros(length(h),1);
max_err_RK4 = zeros(length(h),1);
for i = 1 : length(h)
    ti = (t_initial : h(i) : t_final)';
    y_anal = y_analytical(ti);
    [~,y_EE] = EulerExplicitMethod(f,t_initial,t_final,y0,h(i));
    max_err_EE(i) = max(abs(y_EE - y_anal));
    [~,y_EI] = EulerImplicitMethod(f,t_initial,t_final,y0,h(i),epsilon);
    max_err_EI(i) = max(abs(y_EI - y_anal));
    [~,y_RK2] = RK2Method(f,t_initial,t_final,y0,h(i));
    max_err_RK2(i) = max(abs(y_RK2 - y_anal));
    [~,y_RK3] = RK3Method(f,t_initial,t_final,y0,h(i));
    max_err_RK3(i) = max(abs(y_RK3 - y_anal));
    [~,y_RK4] = RK4Method(f,t_initial,t_final,y0,h(i));
    max_err_RK4(i) = max(abs(y_RK4 - y_anal));
end
figure;
loglog(h, max_err_EE,  'r--', 'LineWidth', 0.5); hold on;
loglog(h, max_err_EI,  'g-' , 'LineWidth', 0.5); hold on;
loglog(h, max_err_RK2, 'b--', 'LineWidth', 0.5); hold on;
loglog(h, max_err_RK3, 'm--', 'LineWidth', 0.5); hold on;
loglog(h, max_err_RK4, 'k--', 'LineWidth', 0.5);
grid on;
xlabel('Step Size (h)');
ylabel('Maximum Absolute Error');
title('Methods'' Order of Convergence O(h^k)');
legend('Euler''s Explicit Method', ...
    'Euler''s Implicit Method','R-K 2nd Order Method', ...
    'R-K 3rd Order Method (Rev)', 'R-K 4th Order Method', ...
    'Location', 'Best');
% End of script.