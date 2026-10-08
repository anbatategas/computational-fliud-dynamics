clear
clc

% --- CFD Exercise 5.1 ---
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
Nx = 40;
% Equally spaced time    grid (Nt):
Nt = 7000;
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
fprintf('Tmax = %d s\n',L);
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
T_IM = zeros(Nt,Nx);
T_CN = zeros(Nt,Nx);
%
% --- Initial & Boundary Conditions ---
% Initial  condition (t=0, n=1):
T_EX(1,:) = 10;
T_IM(1,:) = 10;
T_CN(1,:) = 10;
% Boundary condition (x=0, i=1):
T_EX(2:Nt,1) = 0;
T_IM(2:Nt,1) = 0;
T_CN(2:Nt,1) = 0;
% Boundary condition (x=L, i=Nx):
T_EX(2:Nt,Nx) = 0;
T_IM(2:Nt,Nx) = 0;
T_CN(2:Nt,Nx) = 0;
%
% --------------------------------
% --- Explicit Scheme Solution ---
% --------------------------------
for n = 1 : Nt-1
    for i = 2 : Nx-1
        T_EX(n+1,i) = lamda*T_EX(n,i-1)+...
                    (1-2*lamda)*T_EX(n,i)+...
                            lamda*T_EX(n,i+1);
    end
end
% ---------------------------------------
% --- Simple Implicit Scheme Solution ---
% ---------------------------------------
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
% --------------------------------------
% --- Crank-Nicolson Scheme Solution ---
% --------------------------------------
% Tridiagonal -1,0,+1 diagonals definition:
a_CN = -lamda*ones(Nx-2,1);         % Lower diagonal (-1).
a_CN(1) = 0;
b_CN = 2*(1+lamda)*ones(Nx-2,1);    % Main  diagonal ( 0).
c_CN = -lamda*ones(Nx-2,1);         % Upper diagonal (+1).
c_CN(Nx-2) = 0;
% Main loop:
for n = 1 : Nt-1
    % Previous time (n) RHS (in column form):
    d_CN = zeros(Nx-2,1);
    for i = 2 : Nx-1
        d_CN(i-1) = lamda*T_CN(n,i-1) + ...
                    2*(1-lamda)*T_CN(n,i) + ...
                                lamda*T_CN(n,i+1);
    end
    % Solution of system for new time (n+1):
    T_CN(n+1,2:Nx-1) = ThomasAlgorithmTriInput(a_CN,b_CN,c_CN,d_CN)';
end
%
% -------------------------------------------
% --- GRAPHS FOR TEMPERATURE DISTRIBUTION ---
% -------------------------------------------
figure('Name', '3D Evolution');
surf(t, x, T_CN');
shading interp; grid on; box on; colormap default;
xlabel('t (s)'); ylabel('x (m)'); zlabel('T(t,x)');
title('Temperature Distribution (Crank-Nicolson)');
view([220 30]);
hold on;
mesh(T,X,T_CN', 'FaceColor','none','LineWidth',0.5, ...
                            'EdgeColor','k','EdgeAlpha',0.5);
hold off;
%
% ---------------------------------------
% --- ARITHMETIC DATA IN TABLE FORMAT ---    {{{{{{{{EDIT}}}}}}}}}}}}}}}}}
% ---------------------------------------
mid_node = ceil(Nx/2);
fprintf('-------------------------------------------\n');
fprintf('--- ARITHMETIC RESULTS (Crank-Nicolson) ---\n');
fprintf('Center Node T(t, x=%.2f)\n', x(mid_node));
fprintf('Time(s) \t Temperature\n');
for n = 1:10
    fprintf('%.4f \t \t %.4f\n', t(n), T_CN(n, mid_node));
end
fprintf('... \t \t ...\n');
for n = Nt-9:Nt
    fprintf('%.4f \t \t %.4f\n', t(n), T_CN(n, mid_node));
end
fprintf('-------------------------------------------\n');
%
% End of script.