%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNT (7,7) of length 15 in Fig.6a

function a_cnt_7_7_15()

addpath ../../Models_App/General/

inp = 'cnt_7_7-15_sgp.out';

SE_model_general( inp );

end