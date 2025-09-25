%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function SE_bench_using_BLinf( inp, inp_inf_OR_ifMob )

restoredefaultpath;
rehash pathreset;

AU2KCAL = 627.50947427719404;

if nargin == 1
    ifInf = false;
    ifMob = false;
elseif nargin >= 2
    if ischar( inp_inf_OR_ifMob ) || isstring( inp_inf_OR_ifMob )
        ifInf = true;
        inp_inf = inp_inf_OR_ifMob;
        ifMob = false;
    elseif islogical( inp_inf_OR_ifMob )
        ifInf = false;
        ifMob = inp_inf_OR_ifMob;
    else
        error( 'Invalid argument inp_inf_OR_ifMob!!!' )
    end
end

fprintf( '\n==========================================================\n' );
if ifInf
    fprintf( 'TESING %s using BL-inf of %s ... \n', inp, inp_inf );
else
    fprintf( 'TESING %s ... \n', inp );
end

if ifInf
    [ noneqBndLenGrp0, noneqBndGrp0 ] = getNoneqCCBondLen( inp );

    Ngrp = length(noneqBndLenGrp0);
    noneqBndLenGrp = getNoneqCCBondLen( inp_inf );
    assert( length(noneqBndLenGrp) == Ngrp)
    bndLen_inf = [];
    bndLen0 = [];
    bond_atIx0 = [];
    for j = 1 : Ngrp
        bondlen = mean( noneqBndLenGrp{j} );
        Nbnd0 = length( noneqBndLenGrp0{j} );
        bndLen_inf = [ bndLen_inf; bondlen*ones(Nbnd0,1) ];

        bndLen0 = [ bndLen0; noneqBndLenGrp0{j}' ];
        bond_atIx0 = [ bond_atIx0; noneqBndGrp0{j} ];
    end


    % Sort bndLen_inf and bndLen0 by the lexicographical order of bond_atIx:
    [ bond_atIx0, ix_sort ] = sortrows( bond_atIx0 );
    bndLen0 = bndLen0( ix_sort, : );
    bndLen_inf = bndLen_inf( ix_sort, : );
end


% Predict bndLend using both bndTypExt and ring-based features:
bndLen_pred = predict_BL_fixHMOrev( inp, ifMob );

% Strain energy:
SE_pred = SE_model_pred( inp, 0, ifMob, bndLen_pred ); % SE using pred. bondlen

fprintf( '\nStrain energy (pred. BL) = %.1f kcal/mol\n', SE_pred );


end