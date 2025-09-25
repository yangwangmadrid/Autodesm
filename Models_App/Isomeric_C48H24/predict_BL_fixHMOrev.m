%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function [ bndLen_pred, bndLen_noneq, bond_atIx ] = ...
    predict_BL_fixHMOrev( inp, ifMob )

restoredefaultpath;
rehash pathreset;

POW_HMOrev = 2.70;

if nargin == 1
    ifMob = false;
end

% Obtained from model_EDFT_EHMO_CCBond_numHH():
%[ POW_HMOrev, R0 ] = model_EDFT_EHMO_CCBond_numHH( POW_HMOrev, false );
R0 = mean([1.39662275 1.39662275 1.39662275 1.39662275 1.39662200 1.39662200]);


% Optimized parameters from ../../Models_Train/BondLength_Models/Isomeric_C48H24/
C = [
    -0.148137756599480
    1.332208922130160
    1.330563662692855
    1.330996221227133
    1.327862258414070
    1.335596424431933
    1.330599847432850
    1.335924634074002
    1.336839517306866
    1.340696793227426
    1.331763086411673
    1.341587895270984
    1.341644393415186
    1.338769484043799
    1.335390507218745
    0.001581690634738
    0.000568469645274
    -0.004900478959216
    -0.002639520053261
    -0.004117268262726
    -0.006266561409960
    -0.012229442115753
    -0.005167175620422
    -0.012905924086209
    -0.009830970498260
    -0.011916366575362
    0.000000000000000
    0.000000000000000
    -0.001332782215604
    0.000000000000000
    -0.001142921716205
    0.011442803419505
    0.004517784743381
    0.003460358686730
    0.003834118382835
    -0.001663184236109
    -0.001516066342591
    -0.001817127650116
    0.000000000000000
    -0.002995549143613
    -0.000263607741078
    -0.004774923269219
    -0.005420125432773
    -0.004125034415421
    -0.002949715919132
    0.002305578761833
    -0.000523526424281
    -0.001170840880369
    -0.001299315232777
    -0.001514585820361
    -0.001893697614190
    -0.001589016575325
    -0.001208403171334
    -0.001924899018512
    -0.004635738384436
    -0.004748122045493
    -0.003188935307116
    -0.003385595684461
    -0.004044343443834
    -0.003070482422433
    -0.003621166767606
    -0.005318237146257
    -0.004859306805426
    ];

load( '../../Models_Train/BondLength_Models/Isomeric_C48H24/bndTyp_tags.mat', ...
    'tag_bndTypExt', 'tag_bndTypRB' );

[ data, bond_atIx, lm, twist_bonds ] = ...
    getData_HMOrev( inp, POW_HMOrev, R0, ifMob );
[ bndOrd, bndTyp, bndTypExt, bndTypRB, bndLen ] = unpackStructOut( data );


N = length( bndLen );
fprintf( 'A total of %i data points (CC bonds)\n', N );

logBO = log( bndOrd );

feat_bndTypExt = encode_with_labels( tag_bndTypExt, bndTypExt );
feat_bndTypRB = encode_with_labels( tag_bndTypRB, bndTypRB );
X = [ logBO, feat_bndTypExt, feat_bndTypRB ];

bndLen_pred = X*C;

bndLen_noneq = bndLen;

MAX_NLOOP = 100;
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

end