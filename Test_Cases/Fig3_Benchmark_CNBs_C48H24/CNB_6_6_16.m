%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% CNB (6,6)-16 in Fig.3

function CNB_6_6_16()

addpath ../../Models_App/General/

inp = '6_6-CNB.out';

SE_model_general( inp );

end