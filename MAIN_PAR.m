% =========================================================================
%             J A I   M A H A D E V   G U R U J I   K R U P A
% =========================================================================
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 4 NOVEMBER 2025

% =========================================================================
%                        S E T U P  S C R I P T
% =========================================================================
% CLEAR COMMAND WINDOW AND WORKSPACE
clear;clc;

DIR     = 'PROB_01\';
INPFILE = 'INP_FILE.yaml';
GFILE   = 'GAUSS.g';

% =========================================================================
%                D R I V E R  ( D O  N O T  E D I T )
% =========================================================================

INP_YAML = ReadYaml([DIR INPFILE]); % READ INPUT FILE
disp('INP. FILE READ >> COMPLETE');

% BEAM STIFFNESS VALUES
FLD.AXX = INP_YAML.AXX;
FLD.BXX = INP_YAML.BXX;
FLD.DXX = INP_YAML.DXX;
FLD.SXX = INP_YAML.SXX;

% DISTRIBUTED LOADS
FLD.FX = INP_YAML.FX;
FLD.QZ = INP_YAML.QZ;

% MESH INFORMATION
NEM = INP_YAML.NEM;
P   = INP_YAML.P;
X0  = INP_YAML.X0;
L   = INP_YAML.L;

% GENERATE 1D MESH
[GLX,NOD,NNM,NPE] = MESH1D(X0,L,NEM,P);
disp('1D MESH GENERATION >> COMPLETE');

% READ ESSENTIAL BOUNDARY CONDITIONS
TEMP = INP_YAML.EBCS.NODES;
TEMP = regexprep(TEMP, 'F\(NNM\)', '@(NNM)');
FHANDLE = str2func(TEMP);
EBC_TABLE(:,1) = FHANDLE(NNM)';
EBC_TABLE(:,2) = (cell2mat(INP_YAML.EBCS.DOFS))';
EBC_TABLE(:,3) = (cell2mat(INP_YAML.EBCS.VALS))';

NDF     = 3;                % DEGREES OF FREEDOM PER NODE
NEQ     = NNM*NDF;          % NUMBER OF GLOBAL EQUATIONS
NET     = NPE*NDF;          % EQUATION PER ELEMENT
METHOD  = INP_YAML.METHOD;  % SPATIAL DISCRETIZATION METHOD
NONLIN  = INP_YAML.NONLIN;  % ANALYSIS OPTION
if(NONLIN>0)
    EPSILON = INP_YAML.EPSILON; % CONVERGENCE TOLERANCE
    ITERMAX = INP_YAML.ITERMAX; % NUMBER OF ITERATIONS BEFORE TERMINATION
    BETA    = INP_YAML.BETA;    % RELAXATION PARAMETER
    NLS     = INP_YAML.NLS;     % NUMBER OF LOAD INCREMENTS
    DP      = ...               % ARRAY OF LOAD INCREMENTS 
    cell2mat(INP_YAML.LOAD_INCREMENTS); 
else
    ITERMAX = 1;
    NLS     = 1;
end

% WE WILL USE EVENLY SPACED NODES
MODE = 0;

% --------------------------------------------------
% PRECOMPUTE SFL AND DSFL AT ALL INTEGRATION POINTS
% --------------------------------------------------
if(strcmp(METHOD, 'FEM'))
    disp("WORK IN PROGRESS! EXITING PROGRAM.");
    return;
else
    NGP  = INP_YAML.NGP;
    [INTERFACE,CDOMAIN,HMAT] = ...
        PRECOMP_DMCDM(P,NGP,GFILE);  
end
% ..................................................


% --------------------------------------------------
% PRECOMPUTE SOME VARIABLES NEEDED FOR IMPOSING THE
% ESSENTIAL BOUNDARY CONDITIONS DURING SPARSE
% ASSEMBLY
% --------------------------------------------------
NSPV      = size(EBC_TABLE,1);
ISPV      = EBC_TABLE(:,1:2);
VSPV      = EBC_TABLE(:,3);
NB        = (ISPV(:,1)-1)*NDF + ISPV(:,2);
IS_DIRICHLET     = false(NEQ,1);
IS_DIRICHLET(NB) = true;
VSPV_MAP         = zeros(NEQ,1);
VSPV_MAP(NB)     = VSPV;
% ..................................................


% --------------------------------------------------
% NATURAL BOUNDARY CONDITIONS MAY BE OPTIONALLY 
% PRESENT. THE FOLLOWING LINES ENSURE THAT THE CODE 
% RUNS EVEN IF THE IMPORTED MESH OBJECT DOES NOT 
% HAVE AN NBC_TABLE DEFINED.
% --------------------------------------------------
if exist('NBC_TABLE', 'var') && ~isempty(NBC_TABLE)
    NSSV = size(NBC_TABLE,1);
else
    NSSV = 0;
end
% ..................................................

% INTIALIZE SOLUTION ARRAYS...
GCU = sparse(NEQ,1); % CURRENT ITERATION SOLUTION
GPU = sparse(NEQ,1); % PREVIOUS ITERATION SOLUTION

disp('PRECOMPUTATION >> COMPLETE');
disp('INITIATING SOLVER...');
if(strcmp(METHOD, 'FEM'))
    disp('METHOD: FEM');
else
    disp('METHOD: DMCDM');
end
disp(['NONLIN: ', num2str(NONLIN)]);
disp('CLOCK -> ON');
tic;

LOAD = 0;

for NL = 1:NLS
    % INITIALIZE DATA FOR CURRENT LOAD STEP
    CONVG = false; % CONVERGENCE FLAG
    ITER  = 0;     % NUMBER OF ITERATIONS
    if(NONLIN>0)
        LOAD = LOAD + DP(NL);
    else
        LOAD = 1;
    end

    % NONLINEAR ITERATIONS
    while ITER<ITERMAX

        % INCREMENT ITERATION COUNTER
        ITER = ITER + 1;
        
        % APPLY RELAXATION
        if(NONLIN>0)
            WCU = GPU*BETA + (1-BETA)*GCU;
        end
    
        % PREALLOCATE CELL ARRAYS TO HOLD PER-ELEMENT ASSEMBLY DATA
        I_STIFF_CELL = cell(NEM,1); % ROW (FOR GLK)
        J_STIFF_CELL = cell(NEM,1); % COL (FOR GLK)
        V_STIFF_CELL = cell(NEM,1); % VAL (FOR GLK)

        I_FORCE_CELL = cell(NEM,1); % ROW (FOR GLF)
        V_FORCE_CELL = cell(NEM,1); % VAL (FOR GLF)

        % ELEMENT CALCULATIONS AND ASSEMBLY
        for N=1:NEM
            % ELEMENT NODES
            NI = NOD(N,:);
            
            % EXTRACT ELXY FROM GLXY
            ELX = GLX(NI);
    
            % LOCAL-TO-GLOBAL DOF MAP
            S = 1;
            INDXS = zeros(1,NET);
            for ND = 1:NDF
                INDXS(S:NDF:end) = NDF*NI - (NDF - S)*ones(1,NPE);
                S = S + 1;
            end
    
            % EXTRACT ELU FROM GPU
            if(NONLIN>0)
                ELEMSOL = WCU(INDXS);
            else
                ELEMSOL = zeros(NET,1);
            end
           
            % CALCULATE ELEMENT STIFFNESS AND FORCE COEFFICIENTS.
            [ELK, ELF] = ESUB_DMCDM_BEAM(FLD, ELX, HMAT,...
                                         INTERFACE,CDOMAIN,NPE,...
                                         NGP,NONLIN,ELEMSOL,...
                                         LOAD);
            
       
            % STORAGE FOR THIS ELEMENT
            LOCAL_I_STIFF = zeros(NET*NET,1);
            LOCAL_J_STIFF = zeros(NET*NET,1);
            LOCAL_V_STIFF = zeros(NET*NET,1);

            LOCAL_I_FORCE = zeros(NET,1);
            LOCAL_V_FORCE = zeros(NET,1);

            STIFF_COUNT = 1;
            FORCE_COUNT = 1;

            % SPARSE FINITE ELEMENT ASSEMBLY + IMPOSITION OF EBCS
            for I = 1:NET
                II = INDXS(I); % GLOBAL ROW INDEX
                if(IS_DIRICHLET(II))
                    % GLOBAL EQUATION TO BE MODIFIED TO INCLUDE EBC
                    LOCAL_I_FORCE(FORCE_COUNT)   = II;
                    LOCAL_V_FORCE(FORCE_COUNT) = VSPV_MAP(II);
                    FORCE_COUNT = FORCE_COUNT + 1;

                    LOCAL_I_STIFF(STIFF_COUNT) = II;
                    LOCAL_J_STIFF(STIFF_COUNT) = II;
                    LOCAL_V_STIFF(STIFF_COUNT) = 1.0;
                    STIFF_COUNT = STIFF_COUNT + 1;
                else
                    % GLOBAL EQUATION TO BE ASSEMBLED AS IS
                    LOCAL_I_FORCE(FORCE_COUNT) = INDXS(I);
                    LOCAL_V_FORCE(FORCE_COUNT) = ELF(I);
                    FORCE_COUNT = FORCE_COUNT + 1;
                    for J = 1:NET
                        JJ = INDXS(J); % GLOBAL COL INDEX
                        LOCAL_I_STIFF(STIFF_COUNT) = II;
                        LOCAL_J_STIFF(STIFF_COUNT) = JJ;
                        LOCAL_V_STIFF(STIFF_COUNT) = ELK(I,J);
                        STIFF_COUNT = STIFF_COUNT + 1;
                    end
                end      
            end

            % TRIM UNUSED ALLOCATIONS
            LOCAL_I_STIFF(STIFF_COUNT:end) = [];
            LOCAL_J_STIFF(STIFF_COUNT:end) = [];
            LOCAL_V_STIFF(STIFF_COUNT:end) = [];

            LOCAL_I_FORCE(FORCE_COUNT:end) = [];
            LOCAL_V_FORCE(FORCE_COUNT:end) = [];

            % STORE IN CELL
            I_STIFF_CELL{N} = LOCAL_I_STIFF;
            J_STIFF_CELL{N} = LOCAL_J_STIFF;
            V_STIFF_CELL{N} = LOCAL_V_STIFF;

            I_FORCE_CELL{N} = LOCAL_I_FORCE;
            V_FORCE_CELL{N} = LOCAL_V_FORCE;

        end % END OF ELEMENT CALCULATIONS LOOP
    
        % CONCATENATE
        LIST_I_STIFF   = vertcat(I_STIFF_CELL{:});
        LIST_J_STIFF   = vertcat(J_STIFF_CELL{:});
        LIST_VAL_STIFF = vertcat(V_STIFF_CELL{:});
        
        LIST_I_FORCE   = vertcat(I_FORCE_CELL{:});
        LIST_VAL_FORCE = vertcat(V_FORCE_CELL{:});
    
        % FORM THE GLOBAL SPARSE MATRICES
        GLK = sparse(LIST_I_STIFF,LIST_J_STIFF,LIST_VAL_STIFF);
        GLF = sparse(LIST_I_FORCE,1,LIST_VAL_FORCE);

        % IMPOSE NATURAL (FORCE) BOUNDARY CONDITIONS (IF ANY)
        if(NSSV>0)
            GLOB_DOFS = (NBC_TABLE(:,1)-1)*NDF + NBC_TABLE(:,2);
            GLF(GLOB_DOFS) = GLF(GLOB_DOFS) + NBC_TABLE(:,3)*LOAD;
        end
    
        % SOLVE SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. 
        SOLU = GLK\GLF;
    
        switch NONLIN
            case 0 % LINEAR ANALYSIS
                GCU = SOLU;
                disp('LINEAR ANALYSIS COMPLETE.');
                break;
            case 1 % FIXED-POINT ITERATION
                GPU = GCU;
                GCU = SOLU;
            case 2 % NEWTON ITERATION
                GPU = GCU;
                GCU = GCU + SOLU;
        end
    
        if(ITER==1)
            if(NONLIN==2)
                VSPV_MAP = zeros(NEQ,1);
            end
        end
    
        % CHECK IF THE ITERATIONS NEED TO BE TERMINATED.
        if(NONLIN>0)
            DELTAU = GCU(:,1) - GPU(:,1);
            NUMER  = dot(DELTAU,DELTAU);
            DENOM  = dot(GCU,GCU);
            ERROR  = sqrt(NUMER/DENOM);
            if(ERROR<=EPSILON)
                disp(['STEP ',num2str(NL)]);
                disp('STATUS: STEP CONVERGED')
                disp(['ITERATIONS = ',num2str(ITER)]);
                CONVG = true;
                break;
            end
        end
        disp(strcat("ITER ",num2str(ITER)));
    
    end
    if(~CONVG && NONLIN>0)
        disp('STATUS: CODE FAILED TO CONVERGE :(');
        disp('TERMINATING CONTINUATION LOOP')
        break;
    end
    if(NONLIN==0)
        break;
    end
    
end
T1 = toc;
disp('CLOCK -> OFF');
disp(['ELAPSED TIME = ',num2str(T1),' SECONDS']);
disp('ENDING DMSOL PROGRAM.')


