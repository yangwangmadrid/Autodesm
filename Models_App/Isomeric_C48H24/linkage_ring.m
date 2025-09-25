%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function lm_ring = linkage_ring( ringcell )

NR = length( ringcell );

lm_ring = zeros( NR, NR );
for j = 1 : NR
    rg1 = ringcell{j};
    for k = j+1 : NR
        rg2 = ringcell{k};

        if ~isempty( intersect( rg1, rg2 ) )
            lm_ring( j,k ) = 1;
            lm_ring( k,j ) = 1;
        end
    end
end

end