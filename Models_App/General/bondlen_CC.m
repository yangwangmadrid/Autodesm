%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function bondlen = bondlen_CC( coord )

MAX_CC_BONDLEN = 1.75;

N = size( coord, 1 );
bondlen = [];
for j = 1 : N
    for k = j+1 : N
        r = norm( coord(j,:) - coord(k,:) );
        if r < MAX_CC_BONDLEN
            bondlen = [ bondlen; r ];
        end
    end
end

end