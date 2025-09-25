%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Fullerene C60-C1(100) (NAPP=10) in Fig.7

function C60_100()

addpath ../../Models_App/General/

inp = 'C60_0100_sgp.out';

SE_model_general( inp );

end