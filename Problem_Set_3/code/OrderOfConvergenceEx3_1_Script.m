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
max_errors = zeros(length(N_values),1);
dr_values  = zeros(length(N_values),1);
% --- Definition of analytical solution function handle ---
u_analytical = @(r) 2*(1-(r).^2);
% --- Main loop over N_values ---
for k = 1:length(N_values)
    N = N_values(k);
    % --- Definition of computational grid ---
    r  = linspace(0,1,N);
    dr = 1/(N-1);
    dr_values(k) = dr;
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
    % --- Maximum absolute errors ---
    max_errors(k) = max(absErrors);
    %
end
figure('Color', 'w', 'Name', 'Method Order of Convergence');
loglog(max_errors, dr_values, 'ko-', 'LineWidth', 1.5, 'MarkerFaceColor', 'k'); hold on;
% Add Reference Lines for Slope 1 and Slope 2
% We anchor them to the last data point for visual comparison
ref_1st_order = max_errors(end) * (dr_values / dr_values(end)).^1;
ref_2nd_order = max_errors(end) * (dr_values / dr_values(end)).^2;

loglog(dr_values, ref_1st_order, 'b--', 'LineWidth', 1.2, 'DisplayName', '1st Order Slope (Reference)');
loglog(dr_values, ref_2nd_order, 'r:', 'LineWidth', 1.2, 'DisplayName', '2nd Order Slope (Reference)');

% Formatting
ylabel('Grid Spacing Δr (Log Scale)');
xlabel('Max Absolute Error L_\infty (Log Scale)');
title(['Convergence Plot (Measured Slope = ' num2str(slope, '%.2f') ')']);
legend('Numerical Error', 'Slope = 1 (O(\Delta r))', 'Slope = 2 (O(\Delta r^2))', 'Location', 'NorthWest');
grid on; axis tight;
% Ensure axes are truly logarithmic
set(gca, 'XScale', 'log', 'YScale', 'log');
%
% End of script.