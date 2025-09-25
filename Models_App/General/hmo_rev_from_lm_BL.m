%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Revised HMO with beta depending on bond length

function hmoSol = hmo_rev_from_lm_BL( lm, bndLen, POW, R0, q )

if( nargin == 4 )
    q = 0;
end


N = length( lm );
A = zeros( N, N );
iBond = 0;
for j = 1 : N
    for k = j+1 : N
        if lm(j,k) == 0
            continue
        end
        iBond = iBond + 1;
        R = bndLen( iBond );

        %A(j,k) = exp( -Alpha*(R-R0) );
        A(j,k) = (R/R0)^(-POW);
        if ~isreal( A(j,k) )
            A(j,k)
            pause
        end
        A(k,j) = A(j,k);
    end
end

hmoSol = hmo( A, q );

end
