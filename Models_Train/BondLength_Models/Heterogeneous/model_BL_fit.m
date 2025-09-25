%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Predict bond lengths in (hyperthetical) reference PAH isomer:

function model_BL_fit()


%******************** Retrive all data in the dataset ********************
restoredefaultpath;
rehash pathreset;

% Generate it if it does not exist:
data_file = '../../../Datasets/Benzenoids_C6H6_C96H24/AllData_BL.mat';
load( data_file );

N = length( bndLen );
fprintf( 'A total of %i data points (CC bonds)\n', N );
%************************************************************************


logBO = log( bndOrd );
[ feat_bndTypExt, tag_bndTypExt ] = onehot_rows( bndTypExt );
[ feat_bndTypRB, tag_bndTypRB ] = onehot_rows( bndTypRB );
X = [ logBO, feat_bndTypExt, feat_bndTypRB ];
save( 'bndTyp_tags.mat', 'tag_bndTypExt', 'tag_bndTypRB' );

C = X \ bndLen;

fprintf( '\n******************************\n' );
for j = 1 : length(C)        
    fprintf( '    %.15f\n', C(j) );
end                          
fprintf( '******************************\n' );


fprintf( '\n==================================================\n' );
fprintf( 'rho:    %.15f\n', -C(1) );
for j = 2 : 15
    fprintf( 'w_{k,i} | %2i |  %.15f\n', tag_bndTypExt(j-1), C(j) );
end
ringtype_numbers = 1:12;
ringtype_labels = ["1A", "2A", "2B", "2C", "3A", "3B", "3C", ...
                   "4A", "4B", "4C", "5A", "6A"];
% Create a table:
T = table(ringtype_numbers', ringtype_labels', ...
    'VariableNames', {'Number', 'Label'});
for j = 16 : length(C)  
    %fprintf( '%6i:  %.15f\n', tag_bndTypRB(j-15), C(j) );      
    rgtyp1 = floor( tag_bndTypRB(j-15) / 100 );
    rgtyp2 = tag_bndTypRB(j-15) - rgtyp1*100;
    if rgtyp1 == 0
        fprintf( 'u_{k,i} | %s | %.15f\n', T.Label(rgtyp2), C(j) );
    else
        fprintf( 'u_{k,i} | %s/%s |:  %.15f\n', T.Label(rgtyp1), T.Label(rgtyp2), C(j) );
    end
end                          
fprintf( '==================================================\n' );


bndLen_pred = X*C;


RM = corrcoef( bndLen_pred, bndLen );
R2 = RM(2)^2;
RMSE = sqrt(mean((bndLen_pred - bndLen).^2)); 
fprintf( '\nR2 = %.8f\n', R2 );
fprintf('RMSE = %.4f Angstrom\n', RMSE);

figure;      
set( gcf, 'Color', 'w' );
set( gca, 'FontSize', 14 );
hold on 
box on 
axis equal

plot( bndLen_pred, bndLen, 'o', 'MarkerFaceColor', 'r' )
xlabel( 'Pred. C-C bond length (Å)', 'FontSize', 18 )
ylabel( 'DFT C-C bond length (Å)', 'FontSize', 18 )

xmin = min(bndLen_pred);
xmax = max(bndLen_pred);
xmin = xmin - (xmax-xmin)*0.04;
xmax = xmax + (xmax-xmin)*0.04;
ymin = min(bndLen);
ymax = max(bndLen);
ymin = ymin - (ymax-ymin)*0.04;
ymax = ymax + (ymax-ymin)*0.04;
xymin = min([xmin,ymin]);
xymax = max([xmax,ymax]);
xlim( [ xymin, xymax] )
ylim( [ xymin, xymax] )
line( [ xymin, xymax], [ xymin, xymax] )
xticks( yticks )


end
