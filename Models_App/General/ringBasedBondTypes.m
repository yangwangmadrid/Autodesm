%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Get ring-based CC bond types

function bndTypRB = ringBasedBondTypes( ringcell, rgTyp, bond_atIx )

Nbnd = size( bond_atIx, 1 );
NR = length( ringcell );

bndTypRB = zeros( Nbnd, 1 );
for j = 1 : Nbnd
    at1 = bond_atIx(j,1);
    at2 = bond_atIx(j,2);

    % Get the two rings that share the bond at1-at2:
    rgix = [];
    for r = 1 : NR
        if all( ismember( [at1, at2], ringcell{r} ) )
            rgix = [ rgix, r ];
            if length(rgix) == 2
                break
            end
        end
    end
    if length( rgix ) == 1
        bndTypRB(j) = rgTyp( rgix );
    elseif length( rgix ) == 2
        rgTyp_arr = sort( rgTyp( rgix ) );
        bndTypRB(j) = rgTyp_arr(1)*100 + rgTyp_arr(2);
    end
end

end
