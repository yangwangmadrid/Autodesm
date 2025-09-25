%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Predict strain energy using the general model of autodesmotic reations:

function [ SE_np, SE_predBL ] = SE_model_general( inp, ifMob, twist_bonds )
% inp: Gaussian *.out file using B3LYP/6-311G(d)//B3LYP/6-31G(d)

AU2KCAL = 627.50947427719404;

if nargin == 1
    ifMob = false;
    twist_bonds = [];
end
if nargin == 2
    twist_bonds = [];
end

[~, name ] = fileparts( inp );
fprintf( 'Evaluate strain energy for %s\n', name );

EDFT = loadDFTEnergy( inp );
fprintf( 'Actual DFT energy = %.8f\n', EDFT )

% Coordinates:
[ coord_all, elem ] = loadcoordx( inp ); % Including H atoms
coord_C = [];
coord_H = [];
for k = 1 : length(elem)
    if strcmp(elem{k},'C') || strcmp(elem{k},'6')
        coord_C = [ coord_C; coord_all(k,:) ];
    elseif strcmp(elem{k},'H') || strcmp(elem{k},'1')
        coord_H = [ coord_H; coord_all(k,:) ];
    end
end

% E_HMO:
if ~ifMob
    hmosol = hmo_sim( coord_C );
else
    if ~isempty( twist_bonds )
        fprintf( 'User-defined twist_bonds are:\n' );
    end
    hmosol = hmo_sim_mob( coord_C, 0, twist_bonds );
end
EHMO = hmosol.Etot;
fprintf( 'HMO energy = %.10f\n', EHMO );


% Number of C atoms:
NC = size(coord_C,1);
% Number of H atoms:
NH = size(coord_H,1);
if NH > 0
    fprintf( 'C%iH%i:', NC, NH );
else
    fprintf( 'C%i:', NC );
end
assert( NC + NH == size(coord_all,1) )

% CCBondLen:
CCBondLen = bondlen_CC( coord_C );
fprintf( '  %i C-C bonds', length(CCBondLen) );

% distHH:
if NH > 0
    distHH = dist_HH_auto( coord_H );
    nHH = length( find(distHH < 2.05) ); % Exclude long H...H contacts
    fprintf( '  %i H...H close contacts', nHH );
else
    distHH = [];
end
fprintf( '\n' );


% Prediction using the best model:
DEref_model = model_pred( EHMO, CCBondLen, distHH, NC, NH );
% Predict CCBondLen using model:
CCBondLen_pred = predict_BL_HMOrev_SC( inp, ifMob, twist_bonds );

fprintf( '\nUsing HMOrev-model-predicted CC bond lengths:\n' )
%figure
% CCBondLen
% plot( CCBondLen, CCBondLen_pred, 'ro' )
% xlabel( 'Original CC bond length (Angstrom)' )
% ylabel( 'Predicted reference CC bond length (Angstrom)' )

DEref_model_predBL = model_pred( EHMO, CCBondLen_pred, distHH, NC, NH )
load( 'data_refs.mat' );
E_ref = (4*NH-NC)/18.*EDFT_ref1 + (NC-NH)/72.*EDFT_ref2; % in a.u.
EDFT_model = DEref_model/AU2KCAL + E_ref;
EDFT_model_predBL = DEref_model_predBL/AU2KCAL + E_ref;
fprintf( 'Model DFT energy using original CC bonds = %.8f\n', EDFT_model )
fprintf( 'Model DFT energy using model-predicted CC bonds = %.8f\n', EDFT_model_predBL )
% SE_np = ( EDFT - EDFT_model )*AU2KCAL;
% fprintf( '\nStrain energy (nonplanar) = %.1f kcal/mol\n', SE_np )
SE_predBL = ( EDFT - EDFT_model_predBL )*AU2KCAL;
fprintf( '\nStrain energy (pred. BL) = %.1f kcal/mol\n', SE_predBL )
fprintf( '============================================================\n\n' )


end

% Best model:
% model_EDFT_EHMO_CCBond_numHH()
function y_pred = model_pred( EHMO, CCBondLen, distHH, NC, NH )

POW = 0.3496263686;
Req = 1.2067223668;

C = [
    56.303406215855262
    -75.907221988681741
    1.383132975297873
    281.936940193539385
    -72.049728013506765
    4.371419145837001
    ];

N = length( EHMO );

V_CC = sum((CCBondLen - Req).^(-POW));

nHH = zeros(N, 1);
for j = 1:N
    %nHH(j) = length( distHHs{j} );
    % NOTE: cd ../ and run analyze_data_all to see distribution of H...H distances:
    distHH_cut = 2.05; % Angstrom
    nHH(j) = length( find( distHH < distHH_cut ) );
end

%X = [ EHMO, V_CC, nHH, NC, NH, ones(N,1) ];
X = [ EHMO, V_CC, 0, NC, NH, ones(N,1) ];
y_pred = X * C;

end
