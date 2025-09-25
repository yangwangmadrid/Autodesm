%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function [ POW, Req, C ] = model_EDFT_EHMO_CCBond_numHH( ifPlot )

addpath ../

AU2KCAL = 627.50947427719404;

if nargin == 0
    ifPlot = true;
end

data = read_data();
% Unpack data:
[ name_arr, plan, EDFTs, EHMOs, NCs, NHs, ...
    distHHs, CCBondLens, CHBondLens ] = unpackStruct(data);


NMol = length( name_arr );
fprintf( '%i PAH molecules\n', NMol );

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



%========== Model fitting: ==========%
y = EDFTs;

% Energy model for C--C sigma bonds:
% Initial guesses for POW and Req
initial_params = [ -0.4, 1.33 ]; % [POW, Req]

options = optimoptions('fminunc', ...
    'Display', 'iter', ...
    'MaxIterations', 1000, ...
    'Algorithm', 'quasi-newton', ... 
    'OptimalityTolerance', 1e-6);

optimized_params = fminunc(@(params) ...
    compute_rmse(params, y, EHMOs, CCBondLens, distHHs), ...
    initial_params, options);


% Optimized model prediction:
[ y_pred, C ] = model_pred( optimized_params, y, EHMOs, CCBondLens, ...
    distHHs );

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
fprintf('\nModel formula:\n');
fprintf('E = (%f eV)*EHMO + (%f kcal/mol)*V_CC + (%f kcal/mol)*nHH + %f (a.u.)\n', ...
    C(1)*27.21138505, ...
    C(2)*AU2KCAL, C(3)*AU2KCAL, C(4) );
fprintf( '\n******************************\n' );
for j = 1 : length(C)        
    fprintf( '    %.15f\n', C(j) );
end                          
fprintf( '******************************\n' );
fprintf('\nModel evaluation:\n');
fprintf('  R2 = %.6f\n', R2);
fprintf('RMSE = %.4f kcal/mol\n', RMSE*AU2KCAL);

if ~ifPlot
    return
end



% Make plots:
figure;
set( gcf, 'Color', 'w' );
set( gca, 'FontSize', 14 );
hold on
box on
axis equal

plot( y_pred, y, 'ro', 'MarkerFaceColor', 'r' )
xlabel( '\Delta{}E_{ref}^{pred} (a.u.)', 'FontSize', 18 )
ylabel( '\Delta{}E_{ref}^{DFT} (a.u.)', 'FontSize', 18 )

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
    distHHs )
POW = params(1);
Req = params(2);

N = length(y);

V_CC = zeros(N, 1);
for j = 1:N
    bondlen = CCBondLens{j};
    V_CC(j) = sum((bondlen - Req).^(-POW));
    %V_CC(j) = sum((bondlen/Req - 1).^(-POW));
end

nHH = zeros(N, 1);
for j = 1:N
    %nHH(j) = length( distHHs{j} );
    % NOTE: Run analyze_data_all to see distribution of H...H distances:
    distHH_cut = 2.05; % Angstrom
    nHH(j) = length( find( distHHs{j} < distHH_cut ) );
end

X = [ EHMOs, V_CC, nHH, ones(N,1) ];
C = X \ y;
y_pred = X * C;

end

% Helper function to compute RMSE for given parameters
function rmse = compute_rmse( params, y, EHMOs, CCBondLens, distHHs )
y_pred = model_pred( params, y, EHMOs, CCBondLens, distHHs );
rmse = sqrt(mean((y_pred - y).^2));
end