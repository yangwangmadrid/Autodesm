%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Armchair CNBs (n,n) in Fig.4

function armchair_CNBs()

addpath ../../Models_App/General/

n_all = 2 : 2 : 16;
for j = 1 : length( n_all )
    n = n_all(j);
    inp = sprintf( '%i_%i-CNB.out', n, n );
    fprintf( '\n****************************************\n' );
    fprintf( '  Armchair CNB (%i,%i):\n', n, n );
    SE_model_general( inp );
    pause

end


end