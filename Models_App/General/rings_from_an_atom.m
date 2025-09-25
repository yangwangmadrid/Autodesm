%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function rings_from_an_atom( fid, rg, nblist, at )

nb = nblist(at,:);
% Number of neighbors of atom at:
Nnb = length(nb);
if( nb(Nnb) == 0 )
    Nnb = Nnb - 1;
end
% Number of atoms in current rg[]:
nrg = length( rg );
% fprintf( 'Atom: %i  ', at );
% printRing( rg )
% pause

for j = 1 : Nnb
    at_next = nb(j);
    % Check at_next is already in rg[]:
    ix = find( rg == at_next );
    rg_curr = rg;
    
    % Avoid large rings:
    if length(rg_curr) > 7
        continue
    end
    
    if ~isempty(ix)
        % Backward atom:
        if ix == nrg-1
            continue
            % Ring complete:
        elseif ix == 1
            if isMonoCyclic( rg_curr, nblist ) ... % Only print monocyclic
                    && length(rg_curr) <= 7
                printRing( fid, canonicalRingNumbering( rg_curr ) )
            end
            continue
        else
            continue
        end
    else
        % Ring growth:
        %         rg, at_next
        %         nblist( at_next, : )
        %         pause
        rg_curr = [ rg_curr, at_next ];
    end
    
    rings_from_an_atom( fid, rg_curr, nblist, at_next );
end

return

end

function printRing( fid, rg )

%fprintf( fid, 'Ring:  ' );
for j = 1 : length(rg)
    fprintf( fid, ' %i', rg(j) );
end
fprintf( fid, '\n' );

end



function b = isMonoCyclic( rg, nblist )

nrg = length(rg);

for j = 1 : nrg
    at = rg(j);
    if j == 1
        at_1 = rg(nrg);
        at1 = rg(2);
    elseif j == nrg
        at_1 = rg(j-1);
        at1 = rg(1);
    else
        at_1 = rg(j-1);
        at1 = rg(j+1);
    end
    nb = nblist(at,:);
    
    ix = find( nb ~= at_1 & nb ~= at1 );
    rest_nb = nb(ix);
    if ~isempty( find( rest_nb == rg ) )
        b = false;
        return
    end
end

b = true;
return

end


% Canonical numbering of a ring. The rule is as follows:
%   - 1. Keep the first two numbers as small as possible.
%   - 2. Being clockwise or counterclockwise does not matter.
function rg_canon = canonicalRingNumbering( rg )

nrg = length(rg);
ix = find( rg == min(rg) );
if ix == 1
    at_1 = rg(nrg);
else
    at_1 = rg(ix-1);
end
if ix == nrg
    at1 = rg(1);
else
    at1 = rg(ix+1);
end

if at_1 < at1
    isReverse = true;
else
    isReverse = false;
end

if( ~isReverse )
    loop = [ ix+1 : nrg, 1 : ix-1 ];
else
    loop = fliplr( [ ix+1 : nrg, 1 : ix-1 ] );
end
loop = [ ix, loop ];

k = 1;
for j = loop
    rg_canon(k) = rg(j);
    k = k + 1;
end


end