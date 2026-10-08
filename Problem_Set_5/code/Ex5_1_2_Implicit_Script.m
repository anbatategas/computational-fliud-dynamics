clear
clc

% --- CFD Exercise 5.1 ---
%      --- Part 2 ---
%  (Simple Implicit Scheme)
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
Nx = 100;
% Equally spaced time    grid (Nt):
Nt = 100;
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
fprintf('--- Simple Implicit Scheme Stability ---\n');
fprintf('Stability parameter λ = %.4f\n',lamda);
fprintf('Numerical scheme is always stable, independently of λ. \n');
%
% ----------------------------
% --- SOLUTION OF I.B.V.P. ---
% ----------------------------
% --- Initializations ---
% Matrices for node values contain time nodes (n) as rows and
% space nodes (i) as columns.
T_IM = zeros(Nt,Nx);
%
% --- Initial & Boundary Conditions ---
% Initial  condition (t=0, n=1):
T_IM(1,:) = 10;
% Boundary condition (x=0, i=1):
T_IM(2:Nt,1) = 0;
% Boundary condition (x=L, i=Nx):
T_IM(2:Nt,Nx) = 0;
%
% --------------------------------------
% --- Simple Implicit Scheme Solution ---
% --------------------------------------
% Tridiagonal -1,0,+1 diagonals definition:
a_IM = -lamda*ones(Nx-2,1);         % Lower diagonal (-1).
a_IM(1) = 0;                        
b_IM = (1+2*lamda)*ones(Nx-2,1);    % Main  diagonal ( 0).
c_IM = -lamda*ones(Nx-2,1);         % Upper diagonal (+1).
c_IM(Nx-2) = 0;                     

% Main loop:
for n = 1 : Nt-1
    % Previous time (n) temperature (in row form):
    d_IM  = T_IM(n, 2:Nx-1);
    % Solution of system for new time (n+1):
    T_IM(n+1,2:Nx-1) = ThomasAlgorithmTriInput(a_IM,b_IM,c_IM,d_IM)';
    % Note: The implementation of ThomasAlgorithmTriImput.m handles the
    %       fact that d_IM is in row form. The output is in column form,
    %       so we make it in row form with " ' " at the end.
end
%
% -------------------------------------------
% --- GRAPHS FOR TEMPERATURE DISTRIBUTION ---
% -------------------------------------------
str_title1 = 'Temperature Distribution - Simple Implicit Scheme';
str_title2 = 'Temperature Spatial Profile - Simple Implicit Scheme';
str_title3 = 'Temperature Time Profile - Simple Implicit Scheme';
str_sub   = sprintf('(N_x=%d spatial nodes, N_t=%d time nodes, D = %.2f)', Nx, Nt, D);
%
% --------------------------------------------------
% --- Figure 1: Temperature T(t,x) Surface Graph ---
% --------------------------------------------------
figure(1);
surf(t, x, T_IM');
zlim([min(T_IM(:)) 10]); % Dynamic z-limit to show negative oscillations.
shading interp; grid on; box on; colormap default;
xlabel('t (s)'); ylabel('x (m)'); zlabel('T(t,x)');
title({str_title1}, 'FontSize', 12, 'FontWeight', 'bold');
subtitle({str_sub}, 'FontSize', 10, 'FontWeight', 'bold');
view([220 30]);
hold on;
mesh(T, X, T_IM', ...
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
    h = plot(x, T_IM(n,:), ...
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
xlim([0 L]); ylim([-1 12]);
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
    h = plot(t, T_IM(:,i)', ...
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
xlim([0 Tmax]); ylim([0 10]);
%
% ---------------------------------------
% --- ARITHMETIC DATA IN TABLE FORMAT ---
% ---------------------------------------
mid = ceil(Nx/2);
fprintf('-------------------------------------------\n');
fprintf('--- ARITHMETIC RESULTS (Simple Implicit Scheme) ---\n');
fprintf('Central Node T(t, x=%.2f)\n', x(mid));
fprintf('Time(s) \t Temperature\n');
for n = 1:10
    fprintf('%.4f \t \t %.4f\n', t(n), T_IM(n, mid));
end
fprintf('... \t \t ...\n');
for n = Nt-9:Nt
    fprintf('%.4f \t \t %.4f\n', t(n), T_IM(n, mid));
end
fprintf('-------------------------------------------\n');
fprintf('Left Node T(t, x=%.2f)\n', x(1));
fprintf('Time(s) \t Temperature\n');
for n = 1:10
    fprintf('%.4f \t \t %.4f\n', t(n), T_IM(n, 1));
end
fprintf('... \t \t ...\n');
for n = Nt-9:Nt
    fprintf('%.4f \t \t %.4f\n', t(n), T_IM(n, 1));
end
fprintf('-------------------------------------------\n');

%
% End of script.