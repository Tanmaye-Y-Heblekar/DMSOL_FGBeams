% JAI MAHADEV GURUJI KRUPA
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 6 NOVEMBER 2025


clear;clc;

h  = 1;
b  = 1;
L  = 100;
E1 = 30e+06;
E2 = 10e+06;
NU = 0.3;
KS = 5/6;
q0 = 1;
n  = 50;


A0 = b*h;
B0 = b*h^2;
I0 = b*h^3/12;
M = E1/E2;
AXX = E2*A0*(M+n)/(1+n);
BXX = E2*B0*n*(M-1)/(2*(1+n)*(2+n));
DXX = E2*I0*((6+3*n+3*n^2)*M + (8*n+3*n^2+n^3))/(6+11*n+6*n^2+n^3);
SXX = KS*E2*A0/(2*(1+NU))*(M+n)/(1+n);

DXXSTR = AXX*DXX-BXX^2;

disp(strcat("# L/h = ",num2str(L/h)));
disp(strcat("# n = ",num2str(n)));
fprintf('# BEAM STIFFNESS VALUES\n');
fprintf('AXX: %.6e\n', AXX);
fprintf('BXX: %.6e\n', BXX);
fprintf('DXX: %.6e\n', DXX);
fprintf('SXX: %.6e\n', SXX);



u_exact = @(xi) BXX/DXXSTR*q0*L^3/12*(xi-3*xi^2+2*xi^3);
w_exact = @(xi) AXX/DXXSTR*q0*L^4/24*(xi^2-2*xi^3+xi^4) + ...
               1/DXX*q0*L^4/24*(xi-xi^2) + ...
               1/SXX*q0*L^2/2*(xi-xi^2);

disp("Exact vertical displacement @ center");
disp(w_exact(0.5));