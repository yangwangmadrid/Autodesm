%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function HHdist = dist_HH_auto( coord )
% coord: ONLY containing H atoms

MAX_HH_DIST = 3.0;
%MAX_FAR_HH_DIST = 2.3;

N = size( coord, 1 );
HHdist = [];
for j = 1 : N
    for k = j+1 : N
        r = norm( coord(j,:) - coord(k,:) );
        if r < MAX_HH_DIST
            HHdist = [ HHdist; r ];
        end
    end
end

% if min(HHdist) >= MAX_FAR_HH_DIST
%     HHdist = [];
%     return
% end

HHdist = sort( HHdist );

% diff_HHdist = [];
% for j = 2 : length(HHdist)
%     diff_HHdist = [ diff_HHdist; HHdist(j)-HHdist(j-1)];
% end
% 
% ix = find( diff_HHdist == max(diff_HHdist) );
% nHH = ix(1);
% HHdist = HHdist(1:nHH);

end