%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% Determine all rings (3, 4, 5, 6 and 7 membered) from a adhacent matrix

function [ outp, coord, lm ] = getRings( inp )

MAX_CC_BONDLEN = 1.75; % Angstrom

coord = loadcoord_carbon( inp );
lm = linkage( coord, MAX_CC_BONDLEN );

outp = regexprep( inp, '\.[A-z]+$', '\.rings' );
% if exist( outp, 'file' )
%     return
% end

NAt = length(lm);

nblist = lm2nblist( lm );


fid = fopen( outp, 'w' );
for at = 1 : NAt
    rings_from_an_atom( fid, at, nblist, at );
end
fclose( fid );

outptmp = sprintf( '%s.tmp.%.4f', outp, rand() );
cmd = sprintf( 'sort -k1n %s | uniq > %s', outp, outptmp );
system( cmd );
cmd = sprintf( 'mv %s %s', outptmp, outp );
system( cmd );

fprintf( sprintf( 'File %s written\n', outp ) );

end