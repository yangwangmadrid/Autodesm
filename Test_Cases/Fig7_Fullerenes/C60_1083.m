%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Fullerene C60-D2(1083) (NAPP=12) in Fig.7

function C60_1083()

addpath ../../Models_App/General/

inp = 'C60_1083_sgp.out';

SE_model_general( inp );

end