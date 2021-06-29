close all

r0 = 3;
r1 = 2;
n0 = 1;
ns = 1.2;

% A = r1^2/2;
% B = r0^2/2;
% 
% C = r1;
% D = r0;
% 
% E = A-B;
% F = C-D;

A = [(r0^2)/2-(r1^2)/2 r0-r1; r0 1];
B = [ns; n0];

X = linsolve(A,B);

nFunc = @(r) X(1)*r + X(2);

nFunc(r0)
nFunc(r1)