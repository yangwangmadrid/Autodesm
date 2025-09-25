%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function [ E, E_plus_ZPE ] = loadDFTEnergy( inp )
% inp: g16's output

fid = fopen( inp, 'r' );
E = 0;
E_plus_ZPE = 0;
while ~feof(fid)
    tline = fgetl( fid );
    if regexp( tline, 'SCF .*Done:' )
        tline = regexprep( tline, '^.*=', '' );
        tline = regexprep( tline, 'A\.U\..*', '' );
        tline = regexprep( tline, 'a\.u\..*', '' );
        E = str2double( strtrim(tline) );
    elseif regexp( tline, 'Sum of electronic and zero-point Energies=' )
        tline = regexprep( tline, '^.*Energies=', '' );
        E_plus_ZPE = str2double( strtrim(tline) );
    end
end
fclose( fid );

if E == 0 || isnan(E)
    error( 'Unable to find DFT energy in %s', inp );
end

end