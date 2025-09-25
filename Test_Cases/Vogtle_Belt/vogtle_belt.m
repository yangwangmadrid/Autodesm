%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Vogtle belt

function vogtle_belt()

addpath ../../Models_App/General/

inp = '8belt.out';

SE_model_general( inp );

end