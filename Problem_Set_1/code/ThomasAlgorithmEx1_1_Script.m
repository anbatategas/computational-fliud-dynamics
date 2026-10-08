clear
clc

% --- CFD Exercise 1.1 ---
% Solution script (using ThomasAlgorithm.m)
%
% --- Definition of matrix A ---
A = [2 -1 0 0 0;
    1 2 -1 0 0;
    0 2 4 -1 0;
    0 0 0 2 -1;
    0 0 0 1 2]
%
% --- Definition of constant vector d ---
d = [1 2 -1 -2 -1]
% Note: It is ok to define d as a row. The function
% in "ThomasAlgorithm.m" will take care of it.
%
% --- Apply Thomas algorithm and get solution ---
x = ThomasAlgorithm(A,d)
% --- Check validity of solution ---
x1 = A\d'
if(abs(x-x1)<eps)
    disp('Solution via Thomas algorithm and via A\b match!')
end
% End of script.