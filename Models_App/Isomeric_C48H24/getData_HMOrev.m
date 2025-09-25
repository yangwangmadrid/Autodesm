%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Get data for each of the CC bonds in a PAH

function [ data, bond_atIx, lm, twist_bonds ] = ...
    getData_HMOrev( inp, POW, R0, ifMob )


if nargin == 3
    ifMob = false;
end

% NOTE: if inp contains keyword 'mob', then regard it as Mobius system
if ~ifMob
    ifMob = ~isempty( regexpi( inp, 'mob', 'once' ) );
end
if ifMob
    fprintf( '\nThis is a MOBIUS system\n' );
end

% Get rings:
[ rings_file, coord, lm ] = getRings( inp );
ringcell = readRings( rings_file );
NR = length( ringcell );

% Exclude single-ring case (benzene):
if NR < 2
    data = [];
    bond_atIx = [];
    return
end


if ~ifMob
    twist_bonds = [];
    [ hmosol, bndLen0 ] = hmo_rev( coord, POW, R0 );
else
    % Automatically get twist_bonds for DOUBLE-STRANDED CNBs:
    ringcell = getRings_general( inp );
    NR = length( ringcell );
    fprintf( '%i rings in %s\n', NR, inp );
    %ringcell{:}
    lm_ring = linkage_ring( ringcell );
    nblist_ring = neighbors_ring( lm_ring );
    % Interring bond between ring 1 and its 1st neighboring ring:
    interbond = intersect( ringcell{1}, ringcell{nblist_ring{1}(1)} );
    % Find within ring 1 the two atoms linked to interbond:
    twist_bonds = zeros(2,2);
    for at = setdiff( ringcell{1}, interbond )
        for j = 1 : 2
            if lm( at, interbond(j) )
                twist_bonds(j,:) = [ interbond(j), at ];
            end
        end
    end

    [ hmosol, bndLen0 ] = hmo_rev_mob( coord, POW, R0, 0, twist_bonds );
end

N = hmosol.N;
Etot = hmosol.Etot;
Gap = hmosol.Gap;

BO = hmosol.BO;
bond_atIx = BO( :, 1:2 );
bndOrd = BO( :, 3 );
Nbnd = size( BO, 1 );

%nblist = lm2nblist( lm );

CN = sum( lm ); % Coordination numbers
% Determine bond type matrix:
% 1 ==> C2-C2
% 2 ==> C2-C3 or C3-C2
% 3 ==> C3-C3
lm_typ = lm;
for j = 1 : N
    for k = j+1 : N
        if lm(j,k)
            if CN(j) + CN(k) == 5
                lm_typ(j,k) = 2;
            elseif CN(j) == 3 && CN(k) == 3
                lm_typ(j,k) = 3;
            end
            lm_typ(k,j) = lm_typ(j,k);
        end
    end
end

% Get bond types:
% 1 ==> C2-C2
% 2 ==> C2-C3 or C3-C2
% 3 ==> C3-C3
bndTyp = zeros( Nbnd, 1 );
for j = 1 : Nbnd
    bndTyp(j) = lm_typ( BO( j, 1 ), BO( j, 2 ) );
end

nblist = lm2nblist( lm );
% NN bonds for each CC bond:
bonds_NN = NNBonds( nblist, ringcell, BO(:,1:2) );
%bonds_NN{:}
% Ext. bond types:
bndTypExt = extendedBondTypes( bonds_NN, CN );

% Ring types:
rgTyp = ringTypes( ringcell );

% Get ring-based bond types:
bndTypRB = ringBasedBondTypes( ringcell, rgTyp, BO(:,1:2) );


% Get actual DFT bond lengths:
bndLen = [];
for b = 1 : Nbnd
    r = norm( coord(BO(b,1),:) - coord(BO(b,2),:) );
    bndLen = [ bndLen; r ];
end
assert( isequal( bndLen0, bndLen ) )

data = packStruct( bndOrd, bndTyp, bndTypExt, bndTypRB, bndLen );

end




function ring_by_bond = ringsFromBond( ringcell, at1, at2 )

ring_by_bond = [];
for j = 1 : length( ringcell )
    rg = ringcell{j};
    if all( ismember( [at1, at2], rg ) )
        ring_by_bond = [ ring_by_bond; rg ];
    end
end

end


function [ rgBO, rgBT, atIx ] = BOs_in_ring( BO, lm_typ, rg, at1, at2 )

ix1 = find( rg == at1 );
ix2 = find( rg == at2 );

direct = 1; % Forward
ix3 = ix_next( ix1 );
if ix3 == ix2
    ix3 = ix_prev( ix1 );
    direct = 0; % Reverse
end
at3 = rg( ix3 );
%[ at1, at3 ]
BO13 = findBO( BO, at1, at3 );
BT13 = lm_typ( at1, at3 );

if direct
    ix4 = ix_prev( ix2 );
else
    ix4 = ix_next( ix2 );
end
assert( ix4 ~= ix1 )
at4 = rg( ix4 );
%[ at2, at4 ];
BO24 = findBO( BO, at2, at4 );
BT24 = lm_typ( at2, at4 );

if direct
    ix5 = ix_next( ix3 );
else
    ix5 = ix_prev( ix3 );
end
assert( ix5 ~= ix1 )
at5 = rg( ix5 );
%[ at3, at5 ]
BO35 = findBO( BO, at3, at5 );
BT35 = lm_typ( at3, at5 );

if direct
    ix6 = ix_prev( ix4 );
else
    ix6 = ix_next( ix4 );
end
assert( ix6 ~= ix2 )
at6 = rg( ix6 );
%[ at4, at6 ]
BO46 = findBO( BO, at4, at6 );
BT46 = lm_typ( at4, at6 );

%[ at5, at6 ]
BO56 = findBO( BO, at5, at6 );
BT56 = lm_typ( at5, at6 );

rgBO = [ BO13, BO24, BO35, BO46, BO56 ];
rgBT = [ BT13, BT24, BT35, BT46, BT56 ];
atIx = [ at1, at2, at3, at4, at5, at6 ];

end



function ixn = ix_next( ix )
if ix == 6
    ixn = 1;
else
    ixn = ix+1;
end
end

function ixn = ix_prev( ix )
if ix == 1
    ixn = 6;
else
    ixn = ix-1;
end
end


function b = findBO( BO, at1, at2 )
for j = 1 : size( BO, 1 )
    if ismember(at1,BO(j,1:2)) && ismember(at2,BO(j,1:2))
        b = BO(j,3);
        break
    end
end

end
