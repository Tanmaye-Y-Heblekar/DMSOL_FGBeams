% =========================================================================
%          J A I   M A H A D E V   G U R U J I   K R U P A
% =========================================================================
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 5 NOVEMBER 2025

function [ELK,ELF] = ESUB_DMCDM_BEAM(FLD, ELXY, HMAT,...
                                     INTERFACE,CDOMAIN,NPE,...
                                     NGP,NONLIN,ELEMSOL,...
                                     LOAD)
    % INITIALIZE ELEMENT ARRAYS
    NDF = 3; % NUMBER OF DOFS PER NODE
    NET = NPE*NDF; % NUMBER OF ELEMENT EQUATIONS
    ELK = zeros(NET,NET); % ELEMENT STIFFNESS MATRIX
    ELF = zeros(NET,1); % ELEMENT FORCE VECTOR

    % EXTRACT BEAM STIFFNESS VALUES
    AXX = FLD.AXX;
    BXX = FLD.BXX;
    DXX = FLD.DXX;
    SXX = FLD.SXX;

    % EXTRACT LOADS
    FX = FLD.FX; % DISTRIBUTED AXIAL FORCE
    QZ = FLD.QZ; % DISTRIBUTED TRANSVERSE FORCE

    % ELEMENT SOLUTION
    ELU = ELEMSOL(1:3:end);
    ELW = ELEMSOL(2:3:end);
    ELS = ELEMSOL(3:3:end);

    P = NPE - 1; % ELEMENT ORDER
    
    % COMPUTE COEFFICIENTS AT CD INTERFACES
    ZEROS_P_NPE = zeros(P,NPE);
    R11 = ZEROS_P_NPE; R21 = ZEROS_P_NPE; R31 = ZEROS_P_NPE; 
    R12 = ZEROS_P_NPE; R22 = ZEROS_P_NPE; R32 = ZEROS_P_NPE;
    R13 = ZEROS_P_NPE; R23 = ZEROS_P_NPE; R33 = ZEROS_P_NPE;
    
    % FETCH PRECOMPUTED SHAPE FUNCTION DATA AT INTERFACES
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

        % COMPUTE TANGENT COEFFICIENTS
        if(NONLIN>1)
            DU = dot(GDSFL,ELU);
            DS = dot(GDSFL,ELS);
        end
    end

    % DOMAIN INTEGRATION FOR 32 & 33 AND FORCE TERMS
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
    
    % FORM THE ELEMENT STIFFNESS MATRIX
    ELK(1:NDF:end,1:NDF:end) = HMAT*R11;
    ELK(1:NDF:end,2:NDF:end) = HMAT*R12;
    ELK(1:NDF:end,3:NDF:end) = HMAT*R13;
    ELK(2:NDF:end,1:NDF:end) = HMAT*R21;
    ELK(2:NDF:end,2:NDF:end) = HMAT*R22;
    ELK(2:NDF:end,3:NDF:end) = HMAT*R23;
    ELK(3:NDF:end,1:NDF:end) = HMAT*R31;
    ELK(3:NDF:end,2:NDF:end) = HMAT*R32 + ELK32;
    ELK(3:NDF:end,3:NDF:end) = HMAT*R33 + ELK33;

    % FORM THE ELEMENT FORCE VECTOR
    ELF(1:NDF:end) = ELF1;
    ELF(2:NDF:end) = ELF2;
    ELF(3:NDF:end) = ELF3;

end