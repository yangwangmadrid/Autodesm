%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNT (6,6) of length 18 in Fig.6a

function a_cnt_6_6_18()

addpath ../../Models_App/General/

inp = 'cnt_6_6-18_sgp.out';

SE_model_general( inp );

end
