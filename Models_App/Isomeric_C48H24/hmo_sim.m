%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Revised HMO with beta depending on bond length

function [ hmoSol, bndLen ] = hmo_sim( coord, q )

if( nargin == 1 )
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

        A(j,k) = 1;

        A(k,j) = A(j,k);
    end
end

hmoSol = hmo( A, q );

end
