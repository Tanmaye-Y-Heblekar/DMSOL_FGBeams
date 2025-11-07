% JAI MAHADEV GURUJI KRUPA
% AUTHOR: TANMAYE YASHODAN HEBLEKAR

clear;clc;

h = 1;
L = 100;
A = h*h;
I = h*h^3/12;
E = 30e+06;
NU = 0.3;
KS = 5/6;


AXX = E*A;
BXX = 0;
DXX = E*I;
SXX = KS/(2*(1+NU))*AXX;

disp(strcat('# L/h = ',num2str(L/h)));
fprintf('# BEAM STIFFNESS VALUES\n');
fprintf('AXX: %.6e\n', AXX);
fprintf('BXX: %.6e\n', BXX);
fprintf('DXX: %.6e\n', DXX);
fprintf('SXX: %.6e\n', SXX);
