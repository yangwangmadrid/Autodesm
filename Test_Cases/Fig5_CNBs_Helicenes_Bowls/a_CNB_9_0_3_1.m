%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNB (9,0,3)-1 in Fig.5a

function a_CNB_9_0_3_1()

addpath ../../Models_App/General/

inp = '6x1_m1_CNB_sgp.out';

SE_model_general( inp );

end