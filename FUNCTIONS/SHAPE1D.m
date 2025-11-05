% JAI MAHADEV GURUJI KRUPA
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 28 MAY 2025
% ABOUT: 
% ======
% SHAPE FUNCTIONS AND THEIR DERIVATIVES FOR 
% 1D ELEMENT OF SPECIFIED ORDER 'P'.

function [SFL,DSFL] = SHAPE1D(XI,P,XIS)
    % INITIALIZE...
    SFL  = zeros(P+1,1);
    DSFL = zeros(P+1,1);
    for I = 1:P+1
        % DEFINE DUAL VARIABLE AND SET SEED
        DualXI = struct('V',XI,'D',1);
        DualPSI = struct('V',1,'D',0);
        for K = 1:P+1
            if(I==K)
                continue
            end
            DualTEMP1 = dualPlus(DualXI,struct('V',-XIS(K),'D',0));
            DualTEMP2 = struct('V',XIS(I)-XIS(K),'D',0);
            DualTEMP3 = dualBy(DualTEMP1,DualTEMP2);
            DualPSI = dualTimes(DualPSI,DualTEMP3);
            SFL(I)  = DualPSI.V;
            DSFL(I) = DualPSI.D;
        end
    end
end

% FUNCTIONS FOR AUTOMATIC DIFFERENTIATION
function O = dualPlus(A,B)
    O.V = A.V + B.V;
    O.D = A.D + B.D;
end
function O = dualTimes(A,B)
    O.V = A.V*B.V;
    O.D = A.V*B.D + A.D*B.V;
end
function O = dualBy(A,B)
    O.V = A.V/B.V;
    O.D = (B.V*A.D - A.V*B.D)/(B.V*B.V);
end