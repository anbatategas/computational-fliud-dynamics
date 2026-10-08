clear
clc

% --- CFD Exercise 5.1 ---
%      --- Part 1 ---
%    (Explicit Scheme)
% Solution script (using ThomasAlgorithmTriInput.m)
%
% -----------------------------
% --- 1D DIFFUSION EQUATION ---
% -----------------------------
% Equation: dT/dt = D * d2T/dx2 , D = {(t,x): 0<t<T, 0<x<L}
% Initial Condition  (t=0): T(0,x) = f(x) = 10
% Boundary Condition (x=0): T(t,0) = g(t) = 0
% Boundary Condition (x=L): T(t,L) = h(t) = 0
% -----------------------------
%
% --------------------------------
% --- DEFINITION OF PARAMETERS ---
% --------------------------------
% --- Defintion of equation & boundary parameters ---
L    = 1;
Tmax = 2;
D    = 1;
% --- Definition of number of nodes ---
% Equally spaced spatial grid (Nx):
Nx = 23;
% Equally spaced time    grid (Nt):
Nt = 2000;
% --- Definition of computational grid ---
x = linspace(0,L,Nx);
t = linspace(0,Tmax,Nt);
[T,X] = meshgrid(t,x);
% Measures of spatial and time divisions:
dx = L/(Nx-1);
dt = Tmax/(Nt-1);
% Stability parameter definition:
lamda = D*dt/(dx^2);
% Print grid initialization on terminal:
fprintf('--- Grid Definition ---\n');
fprintf('Number of spatial nodes: Nx = %d\n',Nx);
fprintf('Number of time    nodes: Nt = %d\n',Nt);
fprintf('Spatial grid step size : dx = %d\n',dx);
fprintf('Time    grid step size : dt = %d\n',dt);
fprintf('--- Diffusion Constant ---\n');
fprintf('D = %d\n',D);
fprintf('--- Domain Boundaries ---\n');
fprintf('L    = %d m\n',L);
fprintf('Tmax = %d s\n',Tmax);
% Numerical Scheme Stability Check:
fprintf('--- Explicit Scheme Stability ---\n');
if(lamda >= 0.5)
    fprintf('Stability parameter λ = %.4f >= 0.5\n',lamda);
    fprintf('Explicit Scheme is unstable.\n');
else
    fprintf('Stability parameter λ = %.4f < 0.5\n',lamda);
    fprintf('Explicit Scheme is stable.\n');
end
%
% ----------------------------
% --- SOLUTION OF I.B.V.P. ---
% ----------------------------
% --- Initializations ---
% Matrices for node values contain time nodes (n) as rows and
% space nodes (i) as columns.
T_EX = zeros(Nt,Nx);
%
% --- Initial & Boundary Conditions ---
% Initial  condition (t=0, n=1):
T_EX(1,:) = 10;
% Boundary condition (x=0, i=1):
T_EX(2:Nt,1) = 0;
% Boundary condition (x=L, i=Nx):
T_EX(2:Nt,Nx) = 0;
%
% --------------------------------
% --- Explicit Scheme Solution ---
% --------------------------------
for n = 1 : Nt-1
    for i = 2 : Nx-1
        T_EX(n+1,i) = lamda*T_EX(n,i-1)+(1-2*lamda)*T_EX(n,i)+...
                                                lamda*T_EX(n,i+1);
    end
end
%
% -------------------------------------------
% --- GRAPHS FOR TEMPERATURE DISTRIBUTION ---
% -------------------------------------------
str_title1 = 'Temperature Distribution - Explicit Scheme';
str_title2 = 'Temperature Spatial Profile - Explicit Scheme';
str_title3 = 'Temperature Time Profile - Explicit Scheme';
str_sub   = sprintf(['(N_x=%d spatial nodes, N_t=%d time nodes,' ...
                            'D = %.2f, λ = %.4f)'], Nx, Nt, D, lamda);
%
% --------------------------------------------------
% --- Figure 1: Temperature T(t,x) Surface Graph ---
% --------------------------------------------------
figure(1);
surf(t, x, T_EX');
zlim([min(T_EX(:)) 10]); % Dynamic z-limit to show negative oscillations.
shading interp; grid on; box on; colormap default;
xlabel('t (s)'); ylabel('x (m)'); zlabel('T(t,x)');
title({str_title1}, 'FontSize', 12, 'FontWeight', 'bold');
subtitle({str_sub}, 'FontSize', 10, 'FontWeight', 'bold');
view([220 30]);
hold on;
mesh(T, X, T_EX', ...
     'FaceColor', 'none', 'LineWidth', 0.8, ...
                'EdgeColor', 'k', 'EdgeAlpha', 0.5);
hold off;
%
% ---------------------------------------------------
% --- Figure 2: Temperature spatial profile Graph ---
% ---------------------------------------------------
figure(2);
hold on; box on; grid on;
% Setup Line Colors
line_colors = lines(Nt); 
% Initialize arrays to store Legend Info
legend_handles = []; % To store the plot objects for the legend
legend_strings = {}; % To store the text "t = ..."
for n = 1 : Nt
    % Plot the line with its specific color
    h = plot(x, T_EX(n,:), ...
             'LineWidth', 0.1, ...
             'Color', line_colors(n,:));
    legend_handles(end+1) = h;
    legend_strings{end+1} = sprintf('t = %.3f s', t(n));
end
% Create the Legend Box
lgd = legend(legend_handles, legend_strings, 'Location', 'EastOutside');
title(lgd, 'Time (s)');      % Add a title to the legend box
lgd.FontSize = 3;
% Graph Formatting
xlabel('x (m)'); ylabel('Temperature T(t,x)');
title({str_title2}, 'FontSize', 12, 'FontWeight', 'bold');
subtitle({str_sub}, 'FontSize', 10, 'FontWeight', 'bold');
xlim([0 L]); ylim([-1 10]);
%
% ------------------------------------------------
% --- Figure 3: Temperature time profile Graph ---
% ------------------------------------------------
figure(3);
hold on; box on; grid on;
% Setup Line Colors
line_colors = lines(Nx); 
% Initialize arrays to store Legend Info
legend_handles = []; % To store the plot objects for the legend
legend_strings = {}; % To store the text "t = ..."
for i = 1 : Nx
    % Plot the line with its specific color
    h = plot(t, T_EX(:,i)', ...
             'LineWidth', 0.1, ...
             'Color', line_colors(i,:));
    legend_handles(end+1) = h;
    legend_strings{end+1} = sprintf('x = %.3f m', x(i));
end
% Create the Legend Box
lgd = legend(legend_handles, legend_strings, 'Location', 'EastOutside');
title(lgd, 'Position x (m)');      % Add a title to the legend box
lgd.FontSize = 3;
% Graph Formatting
xlabel('t (s)'); ylabel('Temperature T(t,x)');
title({str_title3}, 'FontSize', 12, 'FontWeight', 'bold');
subtitle({str_sub}, 'FontSize', 10, 'FontWeight', 'bold');
xlim([0 Tmax]); ylim([-1 10]);
%
% ---------------------------------------
% --- ARITHMETIC DATA IN TABLE FORMAT ---
% ---------------------------------------
mid = ceil(Nx/2);
fprintf('-------------------------------------------\n');
fprintf('--- ARITHMETIC RESULTS (Explicit Scheme) ---\n');
fprintf('Central Node T(t, x=%.2f)\n', x(mid));
fprintf('Time(s) \t Temperature\n');
for n = 1:10
    fprintf('%.4f \t \t %.4f\n', t(n), T_EX(n, mid));
end
fprintf('... \t \t ...\n');
for n = Nt-9:Nt
    fprintf('%.4f \t \t %.4f\n', t(n), T_EX(n, mid));
end
fprintf('-------------------------------------------\n');
fprintf('Left Node T(t, x=%.2f)\n', x(1));
fprintf('Time(s) \t Temperature\n');
for n = 1:10
    fprintf('%.4f \t \t %.4f\n', t(n), T_EX(n, 1));
end
fprintf('... \t \t ...\n');
for n = Nt-9:Nt
    fprintf('%.4f \t \t %.4f\n', t(n), T_EX(n, 1));
end
fprintf('-------------------------------------------\n');

%
% End of script.