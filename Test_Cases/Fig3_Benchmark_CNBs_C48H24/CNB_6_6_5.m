%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNB (6,6)-5 in Fig.3

function CNB_6_6_5()

addpath ../../Models_App/General/

inp = '2x2_2_1_1_CNB_sgp.out';

SE_model_general( inp );

end