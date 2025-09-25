%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Corannulene in Fig.5d

function d_corannulene()

addpath ../../Models_App/General/

inp = 'corannulene_sgp.out';

SE_model_general( inp );

end