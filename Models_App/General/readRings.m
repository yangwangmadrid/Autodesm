%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only



function ringcell = readRings( rings_file )

fid = fopen( rings_file, 'r' );

ringcell = {};
k = 1;
while true
    % read a line
    tline = fgetl(fid);
    str = strsplit( tline, ' ' );
    
    rg = [];
    for j = 1 : length(str)
        if ~isempty( str{j} )
            rg = [ rg, str2num(str{j}) ];
        end
    end
    ringcell{k} = rg;
    k = k + 1;
    
    % end of file
    if( feof(fid) )
        break;
    end 
end

fclose( fid );
end