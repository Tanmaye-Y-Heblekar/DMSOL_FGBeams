% =========================================================================
%        J A I   J A I   M A H A D E V   G U R U J I   K R U P A
% =========================================================================
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 4 NOVEMBER 2025

function [INTERFACE,CDOMAIN] = PRECOMP_DMCDM(P,NGP,GFILE)
    NPE = P+1; % NODES PER ELEMENT
    % NOTE: BELOW WE ASSUME EVENLY SPACED NODES
    XIS = linspace(-1,1,NPE); % PARENT ELEMENT NODE LOCATIONS
    DXI = 2/P; % SPACING BETWEEN TWO ADJACENT NODES
    XIN = linspace(-1+0.5*DXI,1-0.5*DXI,P); % INTERFACE LOCATIONS

    % COMPUTE SF AND DSF AT THE CONTROL DOMAIN INTERFACES
    % ---------------------------------------------------

    % INTIALIZE
    INTERFACE.SFL_ARRAY  = zeros(NPE,P);
    INTERFACE.DSFL_ARRAY = zeros(NPE,P);
    % LOOP OVER THE INTERFACES
    for I = 1:P
        XI = XIN(I);
        [SFL,DSFL] = SHAPE1D(XI,P,XIS);
        INTERFACE.SFL_ARRAY(:,I)  = SFL;
        INTERFACE.DSFL_ARRAY(:,I) = DSFL;
    end
    
    % COMPUTE SF AND DSF AT THE DOMAIN INTEGRATION POINTS
    % ---------------------------------------------------

    % INITIALIZE
    CDOMAIN.SFL_ARRAY  = zeros(NPE,NGP,NPE);
    CDOMAIN.DSFL_ARRAY = zeros(NPE,NGP,NPE);

    % CONVENTION:
    % CDOMAIN.SFL_ARRAY(:,J,K) REPRESENTS THE SF VECTOR OF THE J-TH
    % GAUSS POINT OF THE K-TH CONTROL DOMAIN.
    % CDOMAIN.DSFL_ARRAY IS DEFINED IN THE SAME MANNER.

    CDOMAIN.GPTS = zeros(NGP,NPE);
    CDOMAIN.GWTS = zeros(NGP,NPE);

    % CONVENTION:
    % CDOMAIN.GPTS(:,K) REPRESENTS THE GAUSS POINT VECTOR OF THE K-TH
    % CONTROL DOMAIN.
    % CDOMAIN.GWTS(:,K) IS DEFINED IN THE SAME MANNER.

    % GET GAUSS POINTS AND WEIGHT
    [GAUSPT,GAUSWT] = GAUSS(GFILE,NGP);

    % LOOP OVER CONTROL DOMAINS
    XIN_TEMP = [-1,XIN,+1]; % APPEND END NODES TO INTERFACES
    for K = 1:NPE
        A = XIN_TEMP(K);
        B = XIN_TEMP(K+1);
        CDOMAIN.GPTS = (B-A)/2*GAUSPT + (A+B)/2; % SCALED GAUSS POINTS
        CDOMAIN.GWTS = (B-A)/2*GAUSWT; % SCALED GAUSS WEIGHTS
        
        % LOOP OVER GAUSS POINTS
        for J = 1:NGP
            XI = CDOMAIN.GPTS(J);
            [SFL,DSFL] = SHAPE1D(XI,P,XIS);
            CDOMAIN.SFL_ARRAY(:,J,K) = SFL;
            CDOMAIN.DSFL_ARRAY(:,J,K) = DSFL;
        end
    end
end