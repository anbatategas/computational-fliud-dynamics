function [t,y] = EulerImplicitMethod(f, t_initial, t_final, y0, h, epsilon)
% Function EulerImplicitMethod.m
%---------------------------------------------------------------
% Euler/Heun's Implicit method for I.V.P. with step h:
% { y'(t) = f(t,y), y0 = y(ti)}, t_initial<=t<=t_final }
% Usage of trapezoidal rule for numerical integration
% with tolerance epsilon. The method is iterative and
% uses Newton-Raphson method (function NRMethod.m).
%      --------------------------------------
% CAUTION : NRMethod.m must be in the same working directory
%           with this m-function file (!)
%      --------------------------------------
% Input : - function handle f for the slope.
%         - initial time parameter t_initial.
%         - final   time parameter t_final.
%         - initial condition y0.
%         - step size h.
%         - termination parameter epsilon.
%           (use of "weighted" infinity norm)
%         (no. of partition intervals: N = (t_final-t_initial)/h)
% Output: - vector y = sol. of I.V.P.        ((N+1) x 1).
%         - vector t = partition node values ((N+1) x 1).
%---------------------------------------------------------------

% --- Input validity ---
% Ensure function handle f is properly defined:
while(true)
    try
        f(rand(1), rand(1));
        break;
    catch
        error('Function handle f must take two scalar inputs (t,y).');
    end
end
% Ensure t_initial, t_final are properly defined:
if(length(t_initial) ~= 1)
    error('t_initial must be a number');
elseif(length(t_final) ~= 1)
    error('t_final must be a number');
end
% Ensure epsilon is positive:
if(epsilon <= 0)
    error('Termination parameter epsilon must be positive.');
end

% --- Euler's Implicit method implementation ---
% --------------------------------------------
% --- Initialization ---
% Time domain partition nodes:
t = (t_initial : h : t_final)';     % t_intial=t(1), t_final=t(N+1)
N = floor((t_final-t_initial)/h);   % No. of partition intervals N.
% Solution vector intialization:
y    = zeros(N+1,1);
y(1) = y0;                          % Initial condition y(1)=y0.
% --- Main loop ---
for i = 1:N
    % Initialization of temporary loop variables:
    t_old = t(i);
    t_new = t(i+1);
    y_old = y(i);
    % --- Usage of Newton-Raphson method ---
    % --------------------------------------------------------
    % Inital guess using Euler's explicit method:
    y_new_initial = y_old + h * f(t_old, y_old);
    % --------------------------------------------------------
    % Newton-Raphson method for y(i+1):
    % Definition of 1x1 system F(x) = 0 from the defining
    % trapezoidal rule of Implicit Euler's method:
    % y_new = y_old + (h/2)*(f(t_old, y_old) + f(t_new,y_new))
    % --------------------------------------------------------
    F = @(y_new) y_new - y_old - (h/2)*(f(t_old, y_old) + f(t_new,y_new));
    % --------------------------------------------------------
    % Value of y(i+1)   (using NRMethod.m, see Ex1.3):
    % --------------------------------------------------------
    y(i+1) = NRMethod(F, y_new_initial, epsilon);
end
end