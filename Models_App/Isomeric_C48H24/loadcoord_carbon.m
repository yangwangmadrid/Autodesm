%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function [ coord_C, elem_C ] = loadcoord_carbon( inp )

[ coord, elem ] = loadcoord( inp );

coord_C = [];
elem_C = {};
k = 1;
for j = 1 : length( elem )
    if strcmpi( elem{j}, 'C' ) || strcmpi( elem{j}, '6' )
        coord_C = [ coord_C; coord(j,:) ];
        elem_C{k} = 'C';
        k = k + 1;
    end
end

end
