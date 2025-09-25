%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function model_EDFT_EHMO_CCBond_numHH()

addpath ../

AU2KCAL = 627.50947427719404;

data = read_data_all();
% Unpack data:
[ name_arr, plan, EDFT, EHMO, NC, NH, ...
    distHH_arr, CCBondLen_arr, CHBondLen_arr ] = unpackStruct(data);


NMol0 = length( name_arr );
fprintf( '%i PAH molecules\n', NMol0 );

CCBondLen_max = 0;
for j = 1 : NMol0
    if CCBondLen_max < max( CCBondLen_arr{j} );
        CCBondLen_max = max( CCBondLen_arr{j} );
    end
end
CCBondLen_min = Inf;
for j = 1 : NMol0
    if CCBondLen_min > min( CCBondLen_arr{j} );
        CCBondLen_min = min( CCBondLen_arr{j} );
    end
end
fprintf( 'Max(C-C bondlen) = %.4f\n', CCBondLen_max );
fprintf( 'Min(C-C bondlen) = %.4f\n', CCBondLen_min );

%========== Check planarity ==========%
PLANARITY_THRES = 0.02;
ix_nonplan = find( abs(plan) > PLANARITY_THRES );
for j = 1 : length( ix_nonplan )
    name = name_arr{ ix_nonplan(j) };
    fprintf( 'PAH %s is not exactly planar: planarity = %.8f\n', ...
        name, plan( ix_nonplan(j) ) );
end
if ~isempty( ix_nonplan )
    error( 'ABORTED' )
end
%========================================


% Find reference PAHs:
ix_ref1 = []; % benzene
ix_ref2 = []; % 4-coronene
for j = 1 : NMol0
    if strcmp( name_arr{j}, 'benzene' )
        ix_ref1 = j;
    elseif strcmp( name_arr{j}, '4-coronene_reopt-nosymm' )
        ix_ref2 = j;
    end
end
if isempty(ix_ref1)
    error( 'benzene not found')
end
if isempty(ix_ref2)
    error( '4-coronene not found')
end
EDFT_ref1 = EDFT( ix_ref1 ); % benzene
EDFT_ref2 = EDFT( ix_ref2 ); % 4-coronene
EHMO_ref1 = EHMO( ix_ref1 ); % benzene
EHMO_ref2 = EHMO( ix_ref2 ); % 4-coronene
CCBondLen_ref1 = CCBondLen_arr{ ix_ref1 }; % benzene
CCBondLen_ref2 = CCBondLen_arr{ ix_ref2 }; % 4-coronene
% Save data of refs.:
save( 'data_refs.mat', 'EDFT_ref1', 'EDFT_ref2', ...
    'EHMO_ref1', 'EHMO_ref2', 'CCBondLen_ref1', 'CCBondLen_ref2' );


% Prep. dataset (excluding 2 refs):
NMol = NMol0 - 2;
ix_set = setdiff( 1:NMol0, [ix_ref1,ix_ref2] );
assert( NMol == length(ix_set) )
names = name_arr( ix_set );
plans = plan( ix_set );
EDFTs = EDFT( ix_set );
EHMOs = EHMO( ix_set );
NCs = NC( ix_set );
NHs = NH( ix_set );
distHHs = distHH_arr( ix_set );
CCBondLens = CCBondLen_arr( ix_set );
CHBondLens = CHBondLen_arr( ix_set );
% DEref:
assert( all( (4*NHs-NCs) > 0 ) )
DErefs = (EDFTs - (4*NHs-NCs)/18.*EDFT_ref1 - (NCs-NHs)/72.*EDFT_ref2) ...
    * AU2KCAL;

%========== Model fitting: ==========%
y = DErefs;

% Energy model for C--C sigma bonds:
% Initial guesses for POW and Req
initial_params = [ -0.4, 1.33 ]; % [POW, Req]

% Optimize parameters to minimize RMSE
options = optimset('Display', 'iter', 'MaxIter', 1000);
optimized_params = fminsearch(@(params) ...
    compute_rmse( params, y, EHMOs, CCBondLens, distHHs, NCs, NHs ), ...
    initial_params, options);

% Optimized model prediction:
[ y_pred, C ] = model_pred( optimized_params, y, EHMOs, CCBondLens, ...
    distHHs, NCs, NHs );

% Evaluation with optimized parameters
RMSE = sqrt(mean((y_pred - y).^2));
RM = corrcoef(y_pred, y); R2 = RM(2)^2;

% Display results:
fprintf('\nA total of %i PAH molecules in training set\n', length(y));
POW = optimized_params(1);
Req = optimized_params(2);
fprintf('\nOptimized Parameters:\n');
fprintf('POW = %.10f\n', POW);
fprintf('Req = %.10f Å\n', Req);
fprintf('\nModel fomula:\n');
fprintf('E = (%f eV)*EHMO + (%f kcal/mol)*V_CC + (%f kcal/mol)*nHH + (%f kcal/mol)*NC + (%f kcal/mol)*NH\n', ...
    C(1)/AU2KCAL*27.21138505, ...
    C(2), C(3), C(4), ...
    C(5) );
fprintf( '\n******************************\n' );
for j = 1 : length(C)        
    fprintf( '    %.15f\n', C(j) );
end                          
fprintf( '******************************\n' );
fprintf('\nModel evaluation:\n');
fprintf('  R2 = %.6f\n', R2);
fprintf('RMSE = %.4f kcal/mol\n', RMSE);

% % Check consistency for 2n-term coeff. and 2m-term coeff:
% C_2n = C(1)*EHMO_ref1/18 + C(2)*sum((CCBondLen_ref1 - Req).^(-POW))/18 ...
%     - C(1)*EHMO_ref2/72 - C(2)*sum((CCBondLen_ref2 - Req).^(-POW))/72
% C_2m = C(1)*EHMO_ref2/72 + C(2)*sum((CCBondLen_ref2 - Req).^(-POW))/72 ...
%     - C(1)*EHMO_ref1*2/9 - C(2)*sum((CCBondLen_ref1 - Req).^(-POW))*2/9
% fprintf( 'C_2n = %.6f kcal/mol  vs.  C(4) = %.6f kcal/mol\n', C_2n, C(4) )
% fprintf( 'C_2m = %.6f kcal/mol  vs.  C(5) = %.6f kcal/mol\n', C_2m, C(5) )



% Make plots:
figure;
set( gcf, 'Color', 'w' );
set( gca, 'FontSize', 14 );
hold on
box on
axis equal

plot( y_pred, y, 'ro', 'MarkerFaceColor', 'r' )
xlabel( '\Delta{}E_{ref}^{pred} (kcal/mol)', 'FontSize', 18 )
ylabel( '\Delta{}E_{ref}^{DFT} (kcal/mol)', 'FontSize', 18 )

xmin = min(y_pred);
xmax = max(y_pred);
xmin = xmin - (xmax-xmin)*0.04;
xmax = xmax + (xmax-xmin)*0.04;
ymin = min(y);
ymax = max(y);
ymin = ymin - (ymax-ymin)*0.04;
ymax = ymax + (ymax-ymin)*0.04;
xymin = min([xmin,ymin]);
xymax = max([xmax,ymax]);
xlim( [ xymin, xymax] )
ylim( [ xymin, xymax] )
line( [ xymin, xymax], [ xymin, xymax] )
xticks( yticks )


end

function [ y_pred, C ] = model_pred( params, y, EHMOs, CCBondLens, ...
    distHHs, NCs, NHs )
POW = params(1);
Req = params(2);

N = length(y);

V_CC = zeros(N, 1);
for j = 1:N
    bondlen = CCBondLens{j};
    V_CC(j) = sum((bondlen - Req).^(-POW));
    %V_CC(j) = sum((bondlen./Req - 1).^(-POW));
end

nHH = zeros(N, 1);
for j = 1:N
    %nHH(j) = length( distHHs{j} );
    % NOTE: cd ../ and run analyze_data_all to see distribution of H...H distances:
    distHH_cut = 2.05; % Angstrom
    nHH(j) = length( find( distHHs{j} < distHH_cut ) );
end

X = [ EHMOs, V_CC, nHH, NCs, NHs, ones(N,1) ];
C = X \ y;
y_pred = X * C;

end

% Helper function to compute RMSE for given parameters
function rmse = compute_rmse( params, y, EHMOs, CCBondLens, distHHs, NCs, NHs )
y_pred = model_pred( params, y, EHMOs, CCBondLens, distHHs, NCs, NHs );
rmse = sqrt(mean((y_pred - y).^2));
end