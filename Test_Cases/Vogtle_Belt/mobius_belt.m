%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Mobius Vogtle belt

function mobius_belt()

addpath ../../Models_App/General/

inp = '8mobius.out';

SE_model_general( inp, true, [ 1 2; 71 70 ] );

end