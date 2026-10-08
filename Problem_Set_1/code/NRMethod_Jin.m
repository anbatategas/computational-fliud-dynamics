function [x_final, iter] = NRMethod_Jin(F, J, x_initial, epsilon)
% Function NRMethod_Jin.m
%---------------------------------------------------------------
% Newton-Raphson's method for non-linear system F(x) = 0.
% (--- The Jacobian matrix is given analytically ---)
%
% Input : - function handle F for the system of equations.
%         - function handle J for the Jacobian matrix.
%  (Warning: J must be given in the right 3x3 form by the user!)
%         - inital guess x_initial                (n x 1).
%         - termination parameter epsilon         (1 x 1).
%           (use of "weighted" infinity norm)
% Output: - vector x_final = sol. of (F(x)=0)     (n x 1).
%         - integer iter = no. of iterations      (1 x 1).
%---------------------------------------------------------------

% --- System dimension definition ---
n = length(x_initial);

% --- Input validity ---
% Ensure epsilon is positive:
if(epsilon <= 0)
    error('Termination parameter epsilon must be positive.');
end
% Ensure x_initial is in column form :
x_initial = x_initial(:);
% Ensure F and J are properly sized :
while(true)
    try
        F(rand(n,1));
        break;
    catch
        error('Function handle F and x_initial are incompatible.');
    end
end
while(true)
    try
        J(rand(n,n));
        break;
    catch
        error('Function handle J and x_initial are incompatible.');
    end
end

% --- Newton-Raphson method implementation ---
% --------------------------------------------
% --- Temporary and Input variables initialization ---
% Initializations:
x_temp  = x_initial;
x_final = x_temp;
iter    = 0;
% Termination Ratio ("Weighted" infinty norm):
terminationRatio = 1 + epsilon;
% ---------------------------------------------
% --- Main Iterative Loop ---
while(terminationRatio >= epsilon)
    % --- 1. Calculation of F(x(k)) ---
    F_temp = F(x_temp);
    % Ensure column form (F is a function handle):
    F_temp = F_temp(:);
    % --- 2. Calculation of J(x(k)) ---
    J_temp = J(x_temp);
    % --- 3. Update solution and temp. parameters ---
    % Check if J_temp is singular:
    if(det(J_temp) == 0)
        disp('Method terminated because J became singular.');
        break;
    end
    % F_temp is ensured to be in column form, so:
    x_final = x_temp - J_temp\F_temp;
    terminationRatio = norm(x_final-x_temp, "inf")/norm(x_final,"inf");
    iter = iter + 1;
    x_temp = x_final;
end
% ---------------------------------------------
end