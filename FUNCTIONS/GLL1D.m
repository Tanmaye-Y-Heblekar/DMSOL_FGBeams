% JAI MAHADEV GURUJI KRUPA
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 26 MAY 2025
% ABOUT: 
% ======
% FUNCTION THAT GENERATES THE GAUSS-LOBATTO-LEGENDRE (GLL) NODE LOCATIONS
% FOR A SPECIFIED ORDER 'P'

function GLLPTS = GLL1D(P)
    % P = ORDER OF THE (1D) ELEMENT
    syms x
    p     = sym(zeros(P+1,1)); % LEGENDRE POLYNOMIALS
    pp    = sym(zeros(P+1,1)); % FIRST DERIVATIVE OF LEGENDRE POLYNOMIALS
    p(1)  = x^0;
    p(2)  = x;
    pp(1) = x*0;
    pp(2) = x^0;
    n     = 1;
    for k = 3:P+1
        p(k)  = 1/(n+1)*((2*n+1)*x*p(k-1)-n*p(k-2));
        pp(k) = (n+1)*p(k-1) + x*pp(k-1);
        n = n + 1;
    end
    % CREATE ARRAY OF NODAL COORDINATES (GLL POINTS)...
    GLLPTS = double([-1;roots(coeffs(pp(P+1),'All'));1]);
end