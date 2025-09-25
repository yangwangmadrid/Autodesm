%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% predict_BL_HMOrev in a self-consistent manner:

function [ bndLen_pred, bndLen_noneq, bond_atIx ] = ...
    predict_BL_HMOrev_SC( inp, ifMob, twist_bonds )

if nargin == 1
    ifMob = false;
    twist_bonds = [];
end
if nargin == 2
    twist_bonds = [];
end

% Obtained from BondLengths_Predict_HMOrev/model_BL():
%POW_HMOrev = 2.86;
POW_HMOrev = 2.70;
R0 = mean([1.39662275 1.39662275 1.39662275 1.39662275 1.39662200 1.39662200]);

% C = [ % POW_HMOrev = 2.86
%     -0.140872310319592
%     1.334566207213816
%     1.333126713267456
%     1.332191728902011
%     1.331550441644707
%     1.338218817316233
%     1.333920407362727
%     1.338279448676215
%     1.346582823692343
%     1.345527047112812
%     1.349172235706205
%     1.345529710471578
%     1.354738685364971
%     1.354099543927923
%     1.351032421864178
%     0.001898487789044
%     0.000181125472573
%     -0.003631002591679
%     -0.001954536274328
%     -0.002372121971859
%     -0.006802317441800
%     -0.010489555978607
%     -0.005334790032914
%     -0.011267097065833
%     -0.011481013671675
%     -0.011577024843546
%     0.000000000000000
%     0.000000000000000
%     -0.005134305357749
%     -0.004132789245950
%     0.000000000000000
%     -0.001366646429745
%     -0.003188661184273
%     -0.007280240124603
%     -0.008782082095445
%     -0.009904933659485
%     -0.010026980841655
%     -0.011690174880815
%     -0.013150043273552
%     -0.008325343407465
%     -0.002333642960181
%     -0.009074257930992
%     -0.008123336107247
%     -0.005142783851664
%     -0.001411709649043
%     -0.002246932386007
%     -0.011025923703001
%     -0.010131565341386
%     -0.013026277165393
%     -0.012619877734489
%     -0.013267196104175
%     -0.013414822142657
%     -0.015769408927386
%     -0.015070648165299
%     -0.009400570271303
%     -0.014476199714912
%     -0.015799799723936
%     -0.016036987799802
%     -0.014889219289953
%     -0.014012138027214
%     -0.015565374389335
%     -0.014315439885894
%     -0.015719946417972
%     -0.013634461431047
%     -0.014371850562775
%     -0.016626006033520
%     -0.014102672851801
%     -0.015539945738899
%     -0.016325477961239
%     ];

C = [ % POW_HMOrev = 2.70
    -0.144462739934225
    1.319254373298926
    1.317978265412756
    1.317334834379952
    1.315804760245513
    1.322646494018377
    1.318365595258987
    1.322972489501876
    1.344060028379580
    1.343143159231158
    1.332864799546207
    1.329361425967347
    1.338800895300078
    1.338021479421889
    1.335011223454149
    0.015675270192129
    0.013894823925785
    0.010050003828008
    0.011673466560777
    0.011306179626509
    0.006794235506996
    0.003048290574182
    0.008278261799661
    0.002276858296243
    0.002158950548406
    0.001994508920940
    0.000000000000000
    0.000000000000000
    -0.005217527029917
    -0.004213847708763
    0.013932737345410
    0.012537594934337
    0.010702527806279
    0.006651034845354
    0.005074682053505
    0.003928722603236
    0.003804785254527
    0.002148600632742
    0.000671047628550
    0.005344652755440
    -0.002399980633524
    0.004504409885768
    -0.008360443135261
    -0.005326591530245
    0.012434108645154
    0.011560649007113
    0.002886349201680
    0.003681062130367
    0.000855143981081
    0.001184346778624
    0.000496897415995
    0.000371538802602
    -0.001982048585741
    -0.001230041804703
    0.004319202108559
    -0.000739836549234
    -0.002038241499546
    -0.002269418075408
    -0.001106039706674
    -0.000238845912823
    -0.001792887883439
    -0.000559896241101
    -0.001988869639664
    0.000000000000000
    -0.000723997909944
    -0.002949967859207
    -0.000417516976094
    -0.001838688488872
    -0.002616489721736  
];

load( 'bndTyp_tags.mat', 'tag_bndTypExt', 'tag_bndTypRB' );

[ data, bond_atIx, lm ] = getData_HMOrev( inp, POW_HMOrev, R0, ifMob, twist_bonds );

[ bndOrd, bndTyp, bndTypExt, bndTypRB, bndLen ] = unpackStructOut( data );
assert( min( bndOrd ) >= 0 )

N = length( bndLen );
fprintf( 'A total of %i data points (CC bonds)\n', N );


logBO = log( bndOrd );

feat_bndTypExt = encode_with_labels( tag_bndTypExt, bndTypExt );
feat_bndTypRB = encode_with_labels( tag_bndTypRB, bndTypRB );
% %======================================================================
% % ==> 
% % Extend tag_bndTypRB to more patterns bndTypRB by approx.:
% fprintf( '\nWARNING: Using approximate parameter for bndTypRB pattern 709:\n' );
% fprintf( '\nWARNING: Using approximate parameter for bndTypRB pattern 309:\n' );
% fprintf( '\nWARNING: Using approximate parameter for bndTypRB pattern 910:\n' );
% fprintf( '\nWARNING: Using approximate parameter for bndTypRB pattern 707:\n' );
% fprintf( '\nWARNING: Using approximate parameter for bndTypRB pattern 307:\n' );
% % fprintf( '\nWARNING: Using approximate parameter for bndTypRB pattern 1011:\n' );
% tag_bndTypRB_approx = [ tag_bndTypRB; 709; 309; 910; 707; 307 ];
% % 709 ~= mean(509,609)
% C = [ C; mean(C(15+ [ 35, 41 ] )) ];
% % 309 ~= mean(209,409)
% C = [ C; mean(C(15+ [ 22, 31 ] )) ];
% % 910 ~= 810
% C = [ C; mean(C(15+ [ 46 ] )) ];
% % 707 ~= mean(505,606)
% C = [ C; mean(C(15+ [ 32, 39 ] )) ];
% % 307 ~= mean(306)
% C = [ C; mean(C(15+ [ 27 ] )) ];
% % % 1011 ~= mean(811,911)
% % C = [ C; mean(C(15+ [ 47, 50 ] )) ];
% feat_bndTypRB = encode_with_labels( tag_bndTypRB_approx, bndTypRB );
% %======================================================================
X = [ logBO, feat_bndTypExt, feat_bndTypRB ];

bndLen_pred = X*C;

bndLen_noneq = bndLen;


MAX_NLOOP = 500;
TOL = 1e-6;
for iLoop = 1 : MAX_NLOOP
    fprintf( 'Iteration %i ...\n', iLoop );
    % Recalc. HMO BOs using bndLen_pred[]:
    if ~ifMob
        hmosol = hmo_rev_from_lm_BL( lm, bndLen_pred, POW_HMOrev, R0 );
    else
        hmosol = hmo_rev_mob_from_lm_BL( lm, bndLen_pred, POW_HMOrev, R0, ...
            0, twist_bonds );
    end
    % Update X:
    logBO = log( hmosol.BO( :, 3 ) );
    X = [ logBO, feat_bndTypExt, feat_bndTypRB ];
    bndLen_pred_new = X*C;
    diff_bndLen = bndLen_pred_new - bndLen_pred;
    bndLen_pred = bndLen_pred_new; % Update bndLen_pred
    if norm( diff_bndLen ) < TOL && max( abs( diff_bndLen ) ) < TOL
        fprintf( 'Converged (TOL = %.0E)\n', TOL );
        break
    end
end
if iLoop == MAX_NLOOP
    error( 'Bond lengths NOT converged after %i iteration', MAX_NLOOP )
end


end