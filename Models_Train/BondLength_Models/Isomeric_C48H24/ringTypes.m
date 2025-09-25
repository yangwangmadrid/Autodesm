%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Get types of each of the rings by considering its neighoring rings

function rgTyp = ringTypes( ringcell )

NR = length( ringcell );

% Get linkage btw. rings:
lm = zeros( NR, NR );
for j = 1 : NR
    for k = j+1 : NR
        lm(j,k) = length( intersect( ringcell{j}, ringcell{k} ) ) == 2;
        lm(k,j) = lm(j,k);
    end
end
CN = sum( lm );
nblist = lm2nblist( lm, 6 );

rgTyp = zeros( NR, 1 );
for j = 1 : NR
    % Type 1 (1A):
    if CN(j) == 1
        rgTyp(j) = 1;
        continue
    end

    % Type 5 (5A):
    if CN(j) == 5
        rgTyp(j) = 11;
        continue
    end

    % Type 6 (6A):
    if CN(j) == 6
        rgTyp(j) = 12;
        continue
    end

    rg = ringcell{j};
    nb = nblist(j, nblist(j,:) ~= 0);
    rg1 = ringcell{ nb(1) };
    rg2 = ringcell{ nb(2) };
    % Adjacent bonds:
    ab1 = intersect( rg, rg1 );
    ab2 = intersect( rg, rg2 );

    % Type 2--4 (2A, 2B, 2C): 
    if CN(j) == 2
        % Unshared bonds:
        ub = setdiff( rg, [ab1, ab2] );
        % Type 2 (2A):
        if length(ub) == 3
            rgTyp(j) = 2;
            continue
        end
        assert( length(ub) == 2 )
        d_ub = circDist( rg, ub(1), ub(2) );
        % Type 3 (2B):
        if d_ub == 1
            rgTyp(j) = 3;
            continue
        end
        % Type 3 (2C):
        if d_ub == 3
            rgTyp(j) = 4;
            continue
        end
        error( 'IMPOSSIBLE FOR RING TYPE 2--4 (2A, 2B, 2C)!!!' )
    end

    rg3 = ringcell{ nb(3) };
    % Adjacent bonds:
    ab3 = intersect( rg, rg3 );

    % Type 5--7 (3A, 3B, 3C): 
    if CN(j) == 3
        % Unshared bonds:
        ub = setdiff( rg, [ab1, ab2, ab3] );
        % Type 5 (3A):
        if length(ub) == 2
            rgTyp(j) = 5;
            continue
        end
        % Type 6 (3B):
        if length(ub) == 1
            rgTyp(j) = 6;
            continue
        end
        % Type 7 (3C):
        if length(ub) == 0
            rgTyp(j) = 7;
            continue
        end
        error( 'IMPOSSIBLE FOR RING TYPE 5--7 (3A, 3B, 3C)!!!' )
    end
    
    rg4 = ringcell{ nb(4) };
    % Adjacent bonds:
    ab4 = intersect( rg, rg4 );

    % Type 8--10 (4A, 4B, 4C):
    if CN(j) == 4
        % Unshared bonds:
        ub = setdiff( rg, [ab1, ab2, ab3, ab4] );
        % Type 8 (4A):
        if length(ub) == 1
            rgTyp(j) = 8;
            continue
        end

        % Type 9-10 (4B, 4C):
        if length(ub) == 0
            nbrgs = [ rg1; rg2; rg3; rg4 ];
            for r1 = 1 : 4
                nAdj = 0;
                for r2 = setdiff( 1:4, r1 )
                    if ~isempty( intersect( nbrgs(r1,:), nbrgs(r2,:) ) )
                        nAdj = nAdj + 1;
                        if nAdj == 2
                            break
                        end
                    end
                end
                if nAdj == 2
                    break
                end
            end
            if nAdj == 2
                rgTyp(j) = 9;
            else
                rgTyp(j) = 10;
            end
            continue
        end
        error( 'IMPOSSIBLE FOR RING TYPE 8--10 (4A, 4B, 4C)!!!' )
    end

    error( 'IMPOSSIBLE!!!' )
end


end


% Get (smallest) circular distance btw. two atoms in a ring:
function d = circDist( rg, at1, at2 )

ix1 = find( rg == at1 );
ix2 = find( rg == at2 );

d = ix2 - ix1;
if d < 0
    d = -d;
end
if d > length(rg)/2
    d = length(rg) - d;
end

end