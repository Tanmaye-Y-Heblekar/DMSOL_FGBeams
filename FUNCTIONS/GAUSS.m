% JAI MAHADEV GURUJI KRUPA
% AUTHOR: TANMAYE YASHODAN HEBLEKAR
% DATE: 28 MAY 2025
% ABOUT: 
% ======
% READS THE INTEGRATION WEIGHTS AND POINTS FROM THE SPECIFIED .G FILE.

function [GAUSPT,GAUSWT] = GAUSS(GFILE,NGP)
    FILEID  = fopen(GFILE);
    TLINE   = fgetl(FILEID);
    TLINES  = cell(50,1);
    K = 1;
    while ischar(TLINE)
        TLINES{K,1} = TLINE; TLINE = fgetl(FILEID); K = K + 1;
    end
    fclose(FILEID);
    LINE = TLINES{2};
    CONT = textscan(LINE,'%f','Delimiter',',');
    SF_S = CONT{1};
    SF   = SF_S(NGP);
    GAUSPT = zeros(NGP,1);
    GAUSWT = zeros(NGP,1);
    for NI = 1:NGP
        LINE = TLINES{SF};
        CONT = textscan(LINE,'%f %f %f');
        GAUSWT(NI) = CONT{2};
        GAUSPT(NI) = CONT{3};
        SF = SF + 1;
    end
end