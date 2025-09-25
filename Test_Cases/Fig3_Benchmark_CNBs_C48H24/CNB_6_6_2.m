%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNB (6,6)-2 in Fig.3

function CNB_6_6_2()

addpath ../../Models_App/General/

inp = '2x2_1_1_1_1_CNB_sgp.out';

SE_model_general( inp );

end