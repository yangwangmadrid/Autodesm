%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% loadcoord.m
% Load Cartesian coordinates from a file
%
% Last update: Dec-20-2018 by Yang Wang
%              -> Allowing to read Gaussian's output (final geometry)
%
function [ coord, elem, status ] = loadcoordx( inputFileName )
% Parameters:
%   inputFileName <string>: name of the input file which includes coordinates of molecule
%
% Return:
%   coord <N*3 matrix>: coordinates of molecule


fid = fopen( inputFileName, 'r' );
if( fid == -1 )
    fclose(fid);
    error('Cannot open file ''%s''.', inputFileName);
end

status = 0;

% ===== For Gaussian =====
[ coord, elem, status ] = loadgaussian( inputFileName );
if( length(coord) > 0 )
    fclose(fid);
    return;
end
% ========================

k = 1;
isCoordBegin = 0;
while 1
    % read a line
    tline = fgetl(fid);
    % get info from line string
    elemtemp = sscanf(tline, '%s%*f%*f%*f');
    elemtemp = char(elemtemp);
    coordtemp = sscanf(tline, '%*s%f%f%f')';
    
    % check if is a valid line
    if( length(elemtemp) == 0 || length(coordtemp) ~= 3 )
        % end of file
        if( feof(fid) || isCoordBegin )
            break;
        end
        continue;
    end
    isCoordBegin = 1;
    if( strcmp( elemtemp, 'C') == 0 && strcmp( elemtemp, '6') == 0)
        if( strcmp( elemtemp, 'H') == 1 || strcmp( elemtemp, '1') == 1)
%            fprintf('Atom %i: %s\n', k, elemtemp);
        else
%            fclose(fid);
%            error('Atom %i: %s\nThis is not a all-carbon molecule.', k, elemtemp);
        end
    end

    %if( strcmp( elemtemp, 'C') == 1 | strcmp( elemtemp, '6') == 1)
    %if( strcmp( elemtemp, 'H') == 0 & strcmp( elemtemp, '1') == 0 )
        coord(k,:) = coordtemp;        
        elem{k} = elemtemp';
        k = k + 1;
    %end
    
    % end of file
    if( feof(fid) )
        break;
    end 
end

if( isempty( coord ) )
    fclose(fid);
    error('No coordinates available in file ''%s''.', inputFileName);
end



fclose(fid);

end
