%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNT (5,5) of length 2 in Fig.6a

function a_cnt_5_5_2()

addpath ../../Models_App/General/

inp = 'cnt_5_5-2_sgp.out';

SE_model_general( inp );

end
