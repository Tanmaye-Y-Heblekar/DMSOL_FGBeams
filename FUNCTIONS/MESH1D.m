% =========================================================================
%        J A I   J A I   M A H A D E V   G U R U J I   K R U P A
% =========================================================================
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 6 NOVEMBER 2025

function [GLX,NOD,NNM,NPE] = MESH1D(X0,L,NEM,P)
    NPE = P + 1; % NODES PER ELEMENT
    NNM = NEM*NPE - (NEM-1); % NUMBER OF NODES
    GLX = linspace(X0,X0+L,NNM); % NODE COORDINATES
    % GENERATE ELEMENT CONNECTIVITY MATRIX
    NOD=zeros(NEM,NPE);
    K = 1;
    for I=1:NEM
        NOD(I,:) = linspace(K,K+P,NPE);
        K = K+P;
    end
end