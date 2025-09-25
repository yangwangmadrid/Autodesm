%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% [5]Helicene in Fig.5c

function c_5helicene()

addpath ../../Models_App/General/

inp = '5helicene_sgp.out';

SE_model_general( inp );

end