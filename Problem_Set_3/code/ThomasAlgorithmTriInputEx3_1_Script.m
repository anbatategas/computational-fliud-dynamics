clear
clc

% --- CFD Exercise 3.1 ---
%      --- Part 1 ---
% Solution script (using ThomasAlgorithmTriInput.m)
%
% --------------------------------
% --- DEFINITION OF PARAMETERS ---
% --------------------------------
% --- Definition of number of nodes N ---
N = 101
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
% --- Check validity of solution ---
% Definition of auxilary matrices for proper operation
% of "spdiags" function:
aa = zeros(N,1);
cc = zeros(N,1);
for i = 1:N-1
    aa(i) = a(i+1); cc(i+1) = c(i);
end
% System's matrix "reconstruction":
A = spdiags([aa b cc],-1:1,N,N);
A = full(A);
% Solve system by "A\d" method:
u1 = A\d;
if(abs(u-u1)<10^-8)
    disp('Solution via Thomas(a,b,c,d) algorithm and via A\b match!')
end
% Solve system by "Thomas(A,d)" method:
u1T = ThomasAlgorithm(A,d);
if(abs(u-u1T)<10^-8)
    disp('Solution via Thomas(A,d) and via Thomas(a,b,c,d) match')
end
%
% --- Definition of analytical solution function handle ---
u_analytical = @(r) 2*(1-(r).^2);
% --- Definition of analytical solution at partition nodes ---
u_anal = u_analytical(r);
%
% --- Definition of table for data output ---
i = (1:N)';
T = table(i, r(:), u(:), u_anal(:), abs(u(:)-u_anal(:)));
T = renamevars(T,"Var2","r");
T = renamevars(T,"Var3", "u_numerical");
T = renamevars(T,"Var4", "u_analytical");
T = renamevars(T,"Var5", "abs error");
disp(T(1:10,:));
disp('    ........................................................');
disp(T(end-9:end,:));
%
% --- Plot numerical and analytical solution ---
figure;
plot(u_anal, r, 'b-', 'LineWidth', 0.5); hold on;
plot(u, r, 'r.', 'LineWidth', 0.5);
grid on;
ylabel('r/R');
xlabel('u_{x}(r)/V_{average}');
title(['Analytical vs. Numerical Solution (N=', num2str(N), ')']);
legend('Analytical Solution', 'Numerical Solution', 'Location', 'Best');
%
% ------------------------
%      --- Part 2 ---
% --- Calculation of spatially average velocity ---
% --- 1. Analytical calculation ---
Q = 0.16;
D = 1;
Vav_anal = 4*Q/(pi*D^2)
%
% --- 2. Numerical calculation ---
% Definition of integrant:
w = u(:) .* r(:);
% "Numerical" average speed calculation:
Vav_num  = 2*Vav_anal*trapz(r,w)
%
% End of script.