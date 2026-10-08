function [t,y] = RK3Method(f, t_initial, t_final, y0, h)
% Function RK3Method.m
%---------------------------------------------------------------
% Runge-Kutta of 3rd order method for I.V.P. with step h:
% { y'(t) = f(t,y), y0 = y(ti)}, t_initial<=t<=t_final }
%
% Input : - function handle f for the slope.
%         - initial time parameter t_initial.
%         - final   time parameter t_final.
%         - initial condition y0.
%         - step size h
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

% --- Runge-Kutta of 3rd order method implementation ---
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
    k1 = f(t(i),y(i));
    k2 = f(t(i)+h/2,y(i)+(h/2)*k1);
    k3 = f(t(i+1),y(i)-h*k1+2*h*k2);
    y(i+1) = y(i) + (h/6)*(k1+4*k2+k3);
end
end