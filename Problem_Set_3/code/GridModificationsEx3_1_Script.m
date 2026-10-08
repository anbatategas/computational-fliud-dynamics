clear
clc

% --- CFD Exercise 3.1 ---
%      --- Part 3 ---
% Solution script (using ThomasAlgorithmTriInput.m)
%
% --------------------------------
% --- DEFINITION OF PARAMETERS ---
% --------------------------------
% --- Definition of number of nodes N ---
N_values = [10, 20, 50, 100, 200, 500, 1000];
% --- Definition of analytical solution function handle ---
u_analytical = @(r) 2*(1-(r).^2);
% --- Main loop over N_values ---
for k = 1:length(N_values)
    N = N_values(k);
    % --- Definition of computational grid ---
    r  = linspace(0,1,N);
    dr = 1/(N-1);
    % --- Definition of matrix A=diag(a,b,c) and d ---
    % Initializations:
    a = zeros(N,1)  ;
    b = -2*ones(N,1); b(1) =  1; b(N) = 1;
    c = zeros(N,1)  ; c(1) = -1; c(N) = 0;
    d = zeros(N,1)  ;
    % Definitions:
    for i = 2:N-1
        a(i) = 1 - dr/(2*r(i));
        c(i) = 1 + dr/(2*r(i));
        d(i) = -8*(dr)^2;
    end
    %
    % --------------------------
    % --- SOLUTION OF B.V.P. ---
    % --------------------------
    % --- Apply Thomas algorithm and get solution ---
    u = ThomasAlgorithmTriInput(a,b,c,d);
    %
    % --- Definition of analytical solution at partition nodes ---
    u_anal = u_analytical(r);
    %
    % --- Definition of absolute errors ---
    absErrors = abs(u(:)-u_anal(:));
    % --- Definition of table for data output ---
    i = (1:N)';
    T = table(i, r(:), u(:), u_anal(:), absErrors);
    T = renamevars(T,"Var2","r");
    T = renamevars(T,"Var3", "u_numerical");
    T = renamevars(T,"Var4", "u_analytical");
    T = renamevars(T,"absErrors", "abs error");
    sprintf('N=%d ', N)
    disp(T(1:10,:));
    disp('    ........................................................');
    disp(T(end-9:end,:));
    %
    % --- Plot numerical and analytical solution ---
    u_figure = figure();
    plot(u_anal, r, 'b-', 'LineWidth', 0.5); hold on;
    plot(u, r, 'r.', 'LineWidth', 0.5);
    grid on;
    ylabel('r/R');
    xlabel('u_{x}(r)/V_{average}');
    title(['Analytical vs. Numerical Solution (N=', num2str(N), ')']);
    legend('Analytical Solution', 'Numerical Solution', 'Location', 'Best');
    %
    % --- Plots absolute error for each grid node ---
    err_figure = figure();
    plot(absErrors, r, 'r.', 'LineWidth', 0.5);
    grid on;
    ylabel('r/R');
    xlabel('Absolute Errors');
    title(['Absolute Error for each node (N=', num2str(N), ')']);
    legend('Absolute Error', 'Location', 'Best');
    %
    % ------------------------
end
%
% End of script.