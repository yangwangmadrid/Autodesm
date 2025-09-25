%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNT (8,8) of length 2 in Fig.6a

function a_cnt_8_8_2()

addpath ../../Models_App/General/

inp = 'cnt_8_8-2_sgp.out';

SE_model_general( inp );

end