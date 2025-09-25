%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNT (7,4) of length 1 in Fig.6b

function b_cnt_7_4_1()

addpath ../../Models_App/General/

inp = 'cnt_7_4-1_sgp.out';

SE_model_general( inp );

end