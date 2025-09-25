%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Get rings of a nanobelt and order them following the macrocircle

function ringcell = getRings_general( inp )

% Get rings:
rings_file = getRings( inp );
ringcell = readRings( rings_file );
NR = length( ringcell );

end