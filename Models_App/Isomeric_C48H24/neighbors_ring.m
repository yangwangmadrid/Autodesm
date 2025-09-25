%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function nblist_ring = neighbors_ring( lm_ring )

NR = length( lm_ring );

nblist_ring = cell( NR, 1 );
for j = 1 : NR
    nb = find( lm_ring(j,:) ~= 0 );
    nblist_ring{j} = nb;
end

end