%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Fullerene C60-C2(1777) (NAPP=5) in Fig.7

function C60_1777()

addpath ../../Models_App/General/

inp = 'C60_1777_sgp.out';

SE_model_general( inp );

end