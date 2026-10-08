clear
clc

% --- CFD Exercise 4.1 ---
% Solution script (using JacobiMethod.m, GaussSeidelMethod.m)
%
% --------------------------------
% --- DEFINITION OF PARAMETERS ---
% --------------------------------
% --- Definition of number of nodes N per side---
% Equally spaced grid (Nx = Ny = N):
N = 101;
% --- Definition of computational grid ---
x = linspace(0,1,N);
y = linspace(0,1,N);
[X,Y] = meshgrid(x,y);
% dx = dy = h
h = 1/(N-1);
% --- Analytical calculations ---
Q = 0.17;
a = 1; mu = 0.104;
Wav  = Q/(a^2);
dpdz = -(Q*mu)/(0.035144*a^4);
S = (a^(2)*dpdz)/(mu*Wav);
Sh2 = S*h^(2);
fprintf('--- Grid Definition ---\n');
fprintf('Number of nodes: Nx = Ny = N = %d\n',N);
fprintf('--- Source Term ---\n');
fprintf('S = %d\n',S);
%
% --- Definition of matrices A and b (A*w = b) ---
% Initializations:
N_unknowns = (N-2)*(N-2);
b = zeros(N_unknowns,1);
% (Assemble matrix A by keeping indices and values
%  before intialization, for higher computational speed.)
I = zeros(N_unknowns*5,1);  % Row indices.
J = zeros(N_unknowns*5,1);  % Column indices.
V = zeros(N_unknowns*5,1);  % Values at computational nodes.
count = 1;                  % Counter for index vector elements.
% Main loop for definition:
% (scan the grid "row-by-row")
for j = 2 : (N-1)
    for i = 2 : (N-1)
        % --- Index k for current center node (w_k) ---
        % ( k ranges in 1, ... , N-unknowns=(N-2)^2 )
        k = (j-2)*(N-2) + (i-1);
        %
        % Computational shell with 5 points:
        % (w_right + w_left + w_up + w_down - 4*w_center = S*h^2)
        %
        % --- Defintion of matrix A ---
        % Center node (diagonal element):
        I(count) = k; J(count) = k; V(count) = -4;
        count = count + 1;
        % Right neighbour (i+1):
        if(i < (N-1))
            I(count) = k; J(count) = k+1; V(count) = 1;
            count = count + 1;
        end
        % Left  neighbour (i-1):
        if(i > 2)
            I(count) = k; J(count) = k-1; V(count) = 1;
            count = count + 1;
        end
        % Up    neighbour (j+1):
        if(j < (N-1))
            k_up = (j-1)*(N-2) + (i-1);
            I(count) = k; J(count) = k_up; V(count) = 1;
            count = count + 1;
        end
        % Down  neighbour (j-1):
        if(j > 2)
            k_down = (j-3)*(N-2) + (i-1);
            I(count) = k; J(count) = k_down; V(count) = 1;
            count = count + 1;
        end
        %
        % --- Defintion of matrix b ---
        b(k) = Sh2;
    end
end
% Remove non-needed elements from index vectors:
I = I(1:count-1);
J = J(1:count-1);
V = V(1:count-1);
% Define matrix A instanteously:
A = sparse(I,J,V,N_unknowns,N_unknowns);
disp('Matrix created');
%
% --------------------------
% --- SOLUTION OF B.V.P. ---
% --------------------------
% Definition of inital guess and accuracy:
x_initial = 2.0*ones(N_unknowns,1);
epsilon   = 10^-7;
% --- Apply Jacobi method and get solution ---
[w_J, iter_J] = JacobiMethodV(A,b,x_initial, epsilon);
%
% --- Apply Gauss-Seidel method and get solution ---
[w_GS, iter_GS] = GaussSeidelMethodV(A,b,x_initial, epsilon);
%
% --- Check validity of solution ---
% Solve system by "A\d" method:
%w_M = A\b
%
% ------------------------------------
% --- GRID RECONSTRUCTION & GRAPHS ---
% ------------------------------------
% --- 2D Grid reconstruction ---
% Initialization (including zero boundary conditions):
% Jacobi method:
w_J_full = zeros(N,N);
w_J_full(2:N-1, 2:N-1) = reshape(w_J, N-2, N-2)';
%
% Gauss-Seidel method:
w_GS_full = zeros(N,N);
w_GS_full(2:N-1, 2:N-1) = reshape(w_GS, N-2, N-2)';
%
% --- Plot solution ---
% Jacobi method:
str_title_J1 = 'Dimensionless Velocity Profile w(x,y) - Jacobi Method';
str_title_J2 = 'Dimensionless Velocity Contours - Jacobi Method';
str_title_J3 = 'Dimensionless Velocity Profile w(x,y) - Jacobi Method';
str_sub     = sprintf('(N_x=N_y=N=%d nodes, S = %.2f, Accuracy ε = %.1e)', N, S,epsilon);
% --- Maximum & Mean velocities preparation ---
w_mean_val = 1.0;
[max_val,max_vx] = max(w_J_full(:));
[row_max,col_max] = ind2sub(size(w_J_full),max_vx);
x_peak = X(row_max, col_max);
y_peak = Y(row_max, col_max);
% --- Figure 1(A) : Velocity Profile Surface ---
figure('Color','w', 'Position', [100,100,800,600]);
surf(X,Y,w_J_full);
shading interp; box on; grid on; colormap default;
c = colorbar; c.Label.String = 'w(x,y)';
hold on;
mesh(X,Y,w_J_full, 'FaceColor','none', 'LineWidth', 1.5, ...
    'EdgeColor','k','EdgeAlpha', 0.5);
hold off;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('w(x,y)','FontSize',14);
title({str_title_J1;str_sub}, 'FontSize', 12, 'FontWeight', 'bold');
axis([0 1 0 1 0 2]); % Limits to match slide visual
view(-45, 30);       % Adjust view angle
%
% --- Figure 2(A) : Velocity Contours ---
figure('Color','w', 'Position', [150, 150, 800, 600]);
contour3(X, Y, w_J_full, 20, 'LineWidth', 1.2); % 3D Contours
colormap default;
grid on; box on;
c = colorbar; c.Label.String = 'w(x,y)';
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
zlabel('w(x,y)', 'FontSize', 14);
title({str_title_J2; str_sub}, 'FontSize', 12, 'FontWeight', 'bold');
view(-45, 30);
% Label the maximum velocity:
text(x_peak, y_peak, max_val + 0.1, sprintf('Max Velocity: w_{max} \\approx %.3f', max_val), ...
'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k');
%
% --- Figure 3(A) : Mean & Maximum Velocity ---
figure('Color','w', 'Position', [200, 200, 800, 600]);
% Plot Contours for Mean Velocity (w=1):
[C, hC] = contour(X, Y, w_J_full, 'LineWidth', 1.2);
% 2. Identify the contour line for w = 1.0
% We find a point on the grid closest to w=1.0 along the diagonal y=x
diag_indices = 1:N+1:N*N; % Diagonal indices
diag_vals = w_J_full(diag_indices);
% Find index where value crosses 1.0
[~, idx_near_1] = min(abs(diag_vals - 1.0));
x_mean_loc = X(diag_indices(idx_near_1));
y_mean_loc = Y(diag_indices(idx_near_1));

colormap default; colorbar; axis square; grid on;
xlabel('x', 'FontSize', 14); ylabel('y', 'FontSize', 14);
title({str_title_J3; str_sub}, 'FontSize', 12, 'FontWeight', 'bold');
%
% Adding label boxes for Maximum and Mean Velocity points:
% Box 1: Pointing to Maximum Velocity (at the center)
text(x_peak, y_peak, sprintf('\\bf Max Velocity\n w \\approx %.3f', max_val), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'LineWidth', 1);

% Box 2: Pointing to Mean Velocity (w = 1.0)
text(x_mean_loc + 0.1, y_mean_loc - 0.1, sprintf('\\bf Mean Velocity\n w = %.1f', w_mean_val), ...
    'HorizontalAlignment', 'left', 'LineWidth', 1);

% Arrow from the Mean Box to the w=1 contour
line([x_mean_loc + 0.1, x_mean_loc], [y_mean_loc - 0.1, y_mean_loc], ...
    'Color', 'k', 'LineWidth', 1, 'Marker', '.', 'LineStyle', '--');

fprintf('J Maximum dimensionless velocity w: %.4f\n', max_val);
%
% ----------------------------------------------------------------
%
% Gauss-Seidel method:
str_title_GS1 = 'Dimensionless Velocity Profile w(x,y) - Gauss-Seidel Method';
str_title_GS2 = 'Dimensionless Velocity Contours - Gauss-Seidel Method';
str_title_GS3 = 'Dimensionless Velocity Profile w(x,y) - Gauss-Seidel Method';
str_sub     = sprintf('(N_x=N_y=N=%d nodes, S = %.2f, Accuracy ε = %.1e)', N, S,epsilon);
% --- Maximum & Mean velocities preparation ---
w_mean_val = 1.0;
[max_val,max_vx] = max(w_GS_full(:));
[row_max,col_max] = ind2sub(size(w_GS_full),max_vx);
x_peak = X(row_max, col_max);
y_peak = Y(row_max, col_max);
% --- Figure 1(B) : Velocity Profile Surface ---
figure('Color','w', 'Position', [100,100,800,600]);
surf(X,Y,w_GS_full);
shading interp; box on; grid on; colormap default;
c = colorbar; c.Label.String = 'w(x,y)';
hold on;
mesh(X,Y,w_GS_full, 'FaceColor','none', 'LineWidth', 1.5, ...
    'EdgeColor','k','EdgeAlpha', 0.5);
hold off;
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('w(x,y)','FontSize',14);
title({str_title_GS1;str_sub}, 'FontSize', 12, 'FontWeight', 'bold');
axis([0 1 0 1 0 2]); % Limits to match slide visual
view(-45, 30);       % Adjust view angle
%
% --- Figure 2(B) : Velocity Contours ---
figure('Color','w', 'Position', [150, 150, 800, 600]);
contour3(X, Y, w_GS_full, 20, 'LineWidth', 1.2); % 3D Contours
colormap default;
grid on; box on;
c = colorbar; c.Label.String = 'w(x,y)';
xlabel('x', 'FontSize', 14);
ylabel('y', 'FontSize', 14);
zlabel('w(x,y)', 'FontSize', 14);
title({str_title_GS2; str_sub}, 'FontSize', 12, 'FontWeight', 'bold');
view(-45, 30);
% Label the maximum velocity:
text(x_peak, y_peak, max_val + 0.1, sprintf('Max Velocity: w_{max} \\approx %.3f', max_val), ...
'HorizontalAlignment', 'center', 'BackgroundColor', 'w', 'EdgeColor', 'k');
%
% --- Figure 3(B) : Mean & Maximum Velocity ---
figure('Color','w', 'Position', [200, 200, 800, 600]);
% Plot Contours for Mean Velocity (w=1):
[C, hC] = contour(X, Y, w_GS_full, 'LineWidth', 1.2);
% 2. Identify the contour line for w = 1.0
% We find a point on the grid closest to w=1.0 along the diagonal y=x
diag_indices = 1:N+1:N*N; % Diagonal indices
diag_vals = w_GS_full(diag_indices);
% Find index where value crosses 1.0
[~, idx_near_1] = min(abs(diag_vals - 1.0));
x_mean_loc = X(diag_indices(idx_near_1));
y_mean_loc = Y(diag_indices(idx_near_1));

colormap default; colorbar; axis square; grid on;
xlabel('x', 'FontSize', 14); ylabel('y', 'FontSize', 14);
title({str_title_GS3; str_sub}, 'FontSize', 12, 'FontWeight', 'bold');
%
% Adding label boxes for Maximum and Mean Velocity points:
% Box 1: Pointing to Maximum Velocity (at the center)
text(x_peak, y_peak, sprintf('\\bf Max Velocity\n w \\approx %.3f', max_val), ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'LineWidth', 1);

% Box 2: Pointing to Mean Velocity (w = 1.0)
text(x_mean_loc + 0.1, y_mean_loc - 0.1, sprintf('\\bf Mean Velocity\n w = %.1f', w_mean_val), ...
    'HorizontalAlignment', 'left', 'LineWidth', 1);

% Arrow from the Mean Box to the w=1 contour
line([x_mean_loc + 0.1, x_mean_loc], [y_mean_loc - 0.1, y_mean_loc], ...
    'Color', 'k', 'LineWidth', 1, 'Marker', '.', 'LineStyle', '--');

fprintf('GS Maximum dimensionless velocity w: %.4f\n', max_val);
%
%
% ---------------------------------------
% --- ARITHMETIC DATA IN TABLE FORMAT ---
% ---------------------------------------
% 1. Print Convergence Information
fprintf('\n========================================\n');
fprintf('       CONVERGENCE DATA\n');
fprintf('========================================\n');
fprintf('Grid Size (NxN)    : %d x %d\n', N, N);
fprintf('Initial guess      : x_initial = %d * ones(N,1)\n', x_initial(1));
fprintf('Accuracy (epsilon) : %.1e\n', epsilon);
fprintf('----------------------------------------\n');
fprintf('Jacobi Method      : %d iterations\n', iter_J);
fprintf('Gauss-Seidel Method: %d iterations\n', iter_GS);
fprintf('========================================\n');

% 2. Create coordinate vectors for the full grid
[J_idx, I_idx] = meshgrid(1:N, 1:N); % Indices (j corresponds to y, i to x)

% Jacobi method:
% Flatten all matrices into column vectors
Node_Index = (1:N*N)';
i_index    = I_idx(:);     % Grid Index i (x-direction)
j_index    = J_idx(:);     % Grid Index j (y-direction)
x_coord    = X(:);         % Physical coordinate x
y_coord    = Y(:);         % Physical coordinate y
w_value    = w_J_full(:); % Velocity solution w

% Create the Table
ResultsTable = table(Node_Index, i_index, j_index, x_coord, y_coord, w_value, ...
    'VariableNames', {'Node', 'i', 'j', 'x', 'y', 'w_Velocity_J'});

% Display First 10 and Last 10 rows for the Report
fprintf('\n------------------------------------\n');
fprintf('--- SOLUTION DATA TABLE (Jacobi) ---\n');
fprintf('------------------------------------\n');
disp('              --- First 10 Rows ---');
disp(ResultsTable(1:10, :));
disp('     .....................................');
disp('              --- Last 10 Rows ---');
num_rows = height(ResultsTable);
disp(ResultsTable(num_rows-9:num_rows, :));
%
%
% Gauss-Seidel method:
% Flatten all matrices into column vectors
Node_Index = (1:N*N)';
i_index    = I_idx(:);     % Grid Index i (x-direction)
j_index    = J_idx(:);     % Grid Index j (y-direction)
x_coord    = X(:);         % Physical coordinate x
y_coord    = Y(:);         % Physical coordinate y
w_value    = w_GS_full(:); % Velocity solution w

% Create the Table
ResultsTable = table(Node_Index, i_index, j_index, x_coord, y_coord, w_value, ...
    'VariableNames', {'Node', 'i', 'j', 'x', 'y', 'w_Velocity_GS'});

% Display First 10 and Last 10 rows for the Report
fprintf('\n------------------------------------\n');
fprintf('--- SOLUTION DATA TABLE (Gauss-Seidel) ---\n');
fprintf('------------------------------------\n');
disp('              --- First 10 Rows ---');
disp(ResultsTable(1:10, :));
disp('     .....................................');
disp('              --- Last 10 Rows ---');
num_rows = height(ResultsTable);
disp(ResultsTable(num_rows-9:num_rows, :));
%
% End of script.