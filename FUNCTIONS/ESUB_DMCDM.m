% =========================================================================
%          J A I   M A H A D E V   G U R U J I   K R U P A
% =========================================================================
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 5 NOVEMBER 2025

function [ELK,ELF] = ESUB_DMCDM(FLD, ELXY, HMAT,...
                                INTERFACE,CDOMAIN,NPE,...
                                NGP,NONLIN,ELEMSOL,...
                                LOAD)
    ELU = ELEMSOL(1:3:end);
    ELW = ELEMSOL(2:3:end);
    ELS = ELEMSOL(3:3:end);

    P = NPE - 1; % ELEMENT ORDER
    % COMPUTE COEFFICIENTS AT CD INTERFACES
    ZEROS_P_NPE = zeros(P,NPE);
    R11 = ZEROS_P_NPE;
    R12 = ZEROS_P_NPE;
    R13 = ZEROS_P_NPE;
    R21 = ZEROS_P_NPE;
    R22 = ZEROS_P_NPE;
    R23 = ZEROS_P_NPE;
    R31 = ZEROS_P_NPE;
    R32 = ZEROS_P_NPE;
    R33 = ZEROS_P_NPE;

    SFL_ARRAY = INTERFACE.SFL_ARRAY;
    DSFL_ARRAY = INTERFACE.DSFL_ARRAY;
    % LOOP OVER INTERFACES
    for I=1:P
        SFL = SFL_ARRAY(:,I);
        DSFL = DSFL_ARRAY(:,I);
        J = dot(ELXY,DSFL);
        GDSFL = DSFL/J;

        DW = dot(GDSFL,ELW);
        
        R11(I,:) = AXX*GDSFL';
        R12(I,:) = (AXX*DW)/2*GDSFL';
        R13(I,:) = BXX*GDSFL';
        R21(I,:) = AXX*DW*GDSFL';
        R22(I,:) = (SXX+0.5*AXX*DW^2)*GDSFL';
        R23(I,:) = (SXX*SFL'+BXX*DW*GDSFL');
        R31(I,:) = BXX*GDSFL';
        R32(I,:) = (BXX*DW)/2*GDSFL';
        R33(I,:) = DXX*GDSFL';

        if(NONLIN>1)
            DU = dot(GDSFL,ELU);
            DS = dot(GDSFL,ELS);
        end
    end

    % DOMAIN INTEGRATION FOR 32 & 33 TERMS
    ELK32 = zeros(NPE,NPE);
    ELK33 = zeros(NPE,NPE);
    ELF1  = zeros(NPE,1);
    ELF2  = zeros(NPE,1);
    
    % LOOP OVER CONTROL DOMAINS
    for I = 1:NPE
        % GET SCALED INTEGRATION WEIGHTS FOR CURRENT DOMAIN
        GWTS = CDOMAIN.GWTS(:,I);
        
        % INTEGRATION LOOP
        for NG = 1:NGP
            SFL  = CDOMAIN.SFL_ARRAY(:,NG,I);
            DSFL = CDOMAIN.DSFL_ARRAY(:,NG,I);
            J = dot(ELXY,DSFL);
            GDSFL = DSFL/J;
            CNST = J*GWTS(NG);

            ELK32(I,:) = ELK32(I,:) - SXX*GDSFL'*CNST;
            ELK33(I,:) = ELK33(I,:) - SXX*SFL'*CNST;

            ELF1(I) = ELF1(I) + FX*CNST*LOAD;
            ELF2(I) = ELF2(I) + QZ*CNST*LOAD;
        end
    end

    
end