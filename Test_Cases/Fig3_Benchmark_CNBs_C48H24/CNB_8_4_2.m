%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNB (8,4)-2 in Fig.3

function CNB_8_4_2()

addpath ../../Models_App/General/

inp = '2x3_1_1_1_CNB_sgp.out';

SE_model_general( inp );

end