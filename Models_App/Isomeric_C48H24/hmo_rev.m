%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Revised HMO with beta depending on bond length

function [ hmoSol, bndLen ] = hmo_rev( coord, POW, R0, q )

% Alpha = 2;
% %R0 = 1.4; % Angstrom

if( nargin == 1 )
    POW = 0.7;
    %R0 = 1.4; % Angstrom
    R0 = mean([1.39662275 1.39662275 1.39662275 1.39662275 1.39662200 1.39662200]);
    q = 0;
end
if( nargin == 2 )
    %R0 = 1.4; % Angstrom
    R0 = mean([1.39662275 1.39662275 1.39662275 1.39662275 1.39662200 1.39662200]);
    q = 0;
end
if( nargin == 3 )
    q = 0;
end


lm = linkage( coord );
N = length( lm );
A = zeros( N, N );
bndLen = [];
for j = 1 : N
    for k = j+1 : N
        if lm(j,k) == 0
            continue
        end
        R = norm( coord(j,:) - coord(k,:) );
        bndLen = [ bndLen; R ];

        %A(j,k) = exp( -Alpha*(R-R0) );
        A(j,k) = (R/R0)^(-POW);

        A(k,j) = A(j,k);
    end
end

hmoSol = hmo( A, q );

end
