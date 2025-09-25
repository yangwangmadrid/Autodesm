%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function nblist = lm2nblist( lm, NnbMax )

if nargin == 1
    NnbMax = 3;
end

N = length(lm);

nblist = zeros( N, NnbMax ); % Maximum coordination number is 3

for j = 1 : N
    inb = 1;
    for k = 1 : N
        if k == j
            continue
        end

        if( lm(j,k) )
            nblist(j,inb) = k;
            inb = inb + 1;
        end
    end

end