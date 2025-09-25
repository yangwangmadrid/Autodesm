%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% [8]Circulene in Fig.5b

function b_8circulene()

addpath ../../Models_App/General/

inp = '8circulene_sgp.out';

SE_model_general( inp );

end