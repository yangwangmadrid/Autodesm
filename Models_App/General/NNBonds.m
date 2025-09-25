%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Get the four nearest-neighboring (NN) CC bonds surrounding a given CC bond

function bonds_NN = NNBonds( nblist, ringcell, bond_atIx )

Nbnd = size( bond_atIx, 1 );

bonds_NN = cell( Nbnd, 1 );
for j = 1 : Nbnd
    at1 = bond_atIx(j,1);
    at2 = bond_atIx(j,2);

    nb1 = setdiff( nblist( at1, : ), at2 );
    nb1 = nb1( nb1 ~= 0 ); % Remove terminal 0 (if any)
    nb2 = setdiff( nblist( at2, : ), at1 );
    nb2 = nb2( nb2 ~= 0 ); % Remove terminal 0 (if any)

    % Make canonical order btw. at1 and at2:
    if length(nb1) > length(nb2)
        tmp = at1;
        at1 = at2;
        at2 = tmp;
        tmp = nb1;
        nb1 = nb2;
        nb2 = tmp;
    end

    % CASE-I: C2-C2 bonds:
    if length(nb1) == 1 && length(nb2) == 1
        %[at1 at2]
        bond_NN1 = [ at1, nb1 ];
        bond_NN2 = [ at2, nb2 ];
        bonds_NN{j} = [ bond_NN1; bond_NN2 ];
        %bonds_NN{j}
        %pause
        continue
    end

    % CASE-II: C2-C3 bonds:
    if length(nb1) == 1 && length(nb2) == 2 ...
            || length(nb1) == 2 && length(nb2) == 1
        %[at1 at2]
        bond_NN1 = [ at1, nb1 ];
        % Determine the order in nb2:
        if is_co_ring( ringcell, nb1, nb2(1) )
            assert( ~is_co_ring( ringcell, nb1, nb2(2) ) )
            bond_NN2 = [ at2, nb2(1) ];
            bond_NN3 = [ at2, nb2(2) ];
        elseif is_co_ring( ringcell, nb1, nb2(2) )
            bond_NN2 = [ at2, nb2(2) ];
            bond_NN3 = [ at2, nb2(1) ];
        else
            error( 'IMPOSSIBLE FOR C2-C3 BOND!!' );
        end

        bonds_NN{j} = [ bond_NN1; bond_NN2; bond_NN3 ];
        %bonds_NN{j}
        %pause
        continue
    end

    % CASE-III: C3-C3 bonds:
    if length(nb1) == 2 && length(nb2) == 2
        %[at1 at2]
        % Determine order in nb1[]:
        if nblist( nb1(1), 3 ) == 0
            nb1 = nb1([2,1]);
        end
        bond_NN1 = [ at1, nb1(1) ];
        bond_NN2 = [ at1, nb1(2) ];
        % Determine the order in nb2:
        if is_co_ring( ringcell, nb1(1), nb2(1) )
            assert( ~is_co_ring( ringcell, nb1(1), nb2(2) ) )
            assert( ~is_co_ring( ringcell, nb1(2), nb2(1) ) )
            bond_NN3 = [ at2, nb2(1) ];
            bond_NN4 = [ at2, nb2(2) ];
        elseif is_co_ring( ringcell, nb1(1), nb2(2) )
            assert( ~is_co_ring( ringcell, nb1(2), nb2(2) ) )
            bond_NN3 = [ at2, nb2(2) ];
            bond_NN4 = [ at2, nb2(1) ];
        else
            % Switch nb1(1) <--> nb1(2):
            nb1 = nb1([2,1]);
            bond_NN1 = [ at1, nb1(1) ];
            bond_NN2 = [ at1, nb1(2) ];
            % Determine the order in nb2:
            if is_co_ring( ringcell, nb1(1), nb2(1) )
                assert( ~is_co_ring( ringcell, nb1(1), nb2(2) ) )
                assert( ~is_co_ring( ringcell, nb1(2), nb2(1) ) )
                bond_NN3 = [ at2, nb2(1) ];
                bond_NN4 = [ at2, nb2(2) ];
            elseif is_co_ring( ringcell, nb1(1), nb2(2) )
                assert( ~is_co_ring( ringcell, nb1(2), nb2(2) ) )
                bond_NN3 = [ at2, nb2(2) ];
                bond_NN4 = [ at2, nb2(1) ];
            else
                error( 'IMPOSSIBLE FOR C3-C3 BOND!!' );
            end
        end

        bonds_NN{j} = [ bond_NN1; bond_NN2; bond_NN3; bond_NN4 ];
        %[at1 at2], bonds_NN{j}
        %pause
        continue
    end

    error( 'IMPOSSIBLE SCENARIO!!!' )
end

end



function b = is_co_ring( ringcell, at1, at2 )

b = false;
n = 0;
for j = 1 : length( ringcell )
    if ismember( at1, ringcell{j} )
        if ismember( at2, ringcell{j} )
            b = true;
            return
        else
            n = n + 1;
            if n == 3
                return
            end
        end
    end
end

end