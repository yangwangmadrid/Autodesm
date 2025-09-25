%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Best model for strain energy:
% Obtained from model_EDFT_EHMO_CCBond_numHH()

function [ SE, nHH ] = SE_model_pred( inp, numHH_pred, ifMob, CCBondLen )

AU2KCAL = 627.50947427719404;

if nargin == 2
    ifMob = false;
    CCBondLen = [];
elseif nargin == 3
    CCBondLen = [];
end


% Optimized parameters from ../../Models_Train/Energy_Models/Isomeric_C48H24/model_EDFT_EHMO_CCBond_numHH()
POW = -0.4333268764; 
Req = 1.2951595130; % Å
C = [
    0.080315391997920
    0.248592565473450
    0.001425375604011
    -1844.496421707199715
    ];

% E_HMO:
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
% Number of C atoms:
NC = size(coord_C,1);
fprintf( '\nC%i', NC );
% Number of H atoms:
NH = size(coord_H,1);
fprintf( 'H%i\n', NH );
assert( NC + NH == size(coord_all,1) )
assert( NC == 48 && NH == 24 )

% Get E_HMOrev and CC bond lengths:
if ifMob
    % Automatically get twist_bonds for DOUBLE-STRANDED CNBs:
    lm = linkage( coord_C );
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
end

if isempty( CCBondLen )
    if ~ifMob
        [ hmosol, CCBondLen ] = hmo_sim( coord_C );
    else
        [ hmosol, CCBondLen ] = hmo_sim_mob( coord_C, 0, twist_bonds );
    end
else % User's already defined CCBondLen. So don't overload.
    if ~ifMob
        hmosol = hmo_sim( coord_C );
    else
        hmosol = hmo_sim_mob( coord_C, 0, twist_bonds );
    end
end
EHMO = hmosol.Etot;
fprintf( 'E_HMOrev = %.10f\n', EHMO );
fprintf( 'Max(C-C bondlen) = %.4f\n', max(CCBondLen) )
fprintf( 'Min(C-C bondlen) = %.4f\n', min(CCBondLen) )

% V_CC:
V_CC = sum( (CCBondLen - Req).^(-POW) );

% distHH:
distHH = dist_HH_auto( coord_H );
fprintf( 'Max(H...H dist.) = %.4f\n', min(distHH) );

% nHH in actual, strained molecule:
nHH = length( find(distHH < 2.05) ); % Exclude long H...H contacts
fprintf( '%i close H...H contact(s) in actual molecule\n', nHH );

% In hypothetical model strain-free molecule
fprintf( '%i close H...H contact(s) in model molecule\n', numHH_pred );

% Model formula:
% E = (2.090378 eV)*EHMO + (121.303187 kcal/mol)*V_CC ...
%     + (0.660827 kcal/mol)*nHH + -1840.040386 (a.u.)
%
%X = [ EHMO, V_CC, nHH, 1 ];
%X = [ EHMO, V_CC, numHH_pred, 1 ];
X = [ EHMO, V_CC, 0, 1 ];
EDFT_model = X * C;
fprintf( 'Model DFT energy = %.8f\n', EDFT_model )

EDFT = loadDFTEnergy( inp );
fprintf( 'Actual DFT energy = %.8f\n', EDFT )

% % Exclude H...H repulsion energy in actual molecule for it is NOT regarded
% % as strain energy:
% EDFT_corr = EDFT - C(3)*nHH;
% fprintf( 'Actual DFT energy excl. HH repul. = %.8f\n', EDFT_corr )
%SE = ( EDFT_corr - EDFT_model )*AU2KCAL;

SE = ( EDFT - EDFT_model )*AU2KCAL;

end
