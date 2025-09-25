%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% loadgaussian.m
% Read ordinates of molecule from Gaussian output file
function [ coord, elem, status ] = loadgaussian( inputFileName )
% Parameters:
% inputFileName <string>: name of input file
%
% Return:
% coord <N*3 matrix>: coordinates of molecule
% elem <N*1 cell>: element names of molecule

fid=fopen(inputFileName, 'r');
if fid==-1
    fclose(fid);
    error('Canot open file ''%s''.', inputFileName);
end

elem = {};

n = 1; % number of points
k = 1; % number of atoms
coordMarker = 0;
linenum = 1;
isConverged = false;
coord = [];
while 1
    % read a line
    tline=fgetl(fid);    
    linenum = linenum + 1;
    
    if( ~isempty( strfind(tline, 'Normal termination of Gaussian') ) )
        isConverged = true;
    end
    
    % end of file
    if feof(fid)
        break;
    end
    
    % begin the coordinate block
    if( ~isempty( strfind(tline, 'Standard orientation') ) |...
            ~isempty( strfind(tline, 'Input orientation') )|...
            ~isempty( strfind(tline, 'Z-Matrix orientation') ) )
        coordMarker = 1;            
        continue;
    end
    
    % the following four line below the 'Standard orientation' line
    if( coordMarker >= 1 & coordMarker <= 4)
            coordMarker = coordMarker + 1;
            
            % determine g94 or post-g94
            if( coordMarker == 4 )               
                isg04 = isempty( strfind(tline, 'Type') );
            end            
            
            if( coordMarker == 5 )
                k = 1;
            end
            
            continue;
    end
    % now, reach the coordinates line
    if( coordMarker == 5 )
        % get info from line string
        strfield = sscanf(tline, '%f')';
        % check if is a valid line
        if( ~( ( length(strfield) == 5 & isg04 ) | ( length(strfield) == 6 & ~isg04 ) ) )
            if( k == 1 ) % error
                fclose(fid);
                error('Error: the file is corrupted (line %d) ->\n%s',...
                linenum, tline);
                coordMarker = 0;
                continue;
            else % end of coordinate block
                
                % ...
                
                n = n + 1;
                coordMarker = 0;
                continue;
            end
        end

        elementNumber = strfield(2);
	elem{k} = num2str( elementNumber );
%         if( elementNumber ~= 6 )
%             fclose(fid);
%             error('Atom %i: %s\nThis is not a all-carbon molecule.', k, elementNumber);
%         end
        if( isg04)
            coord(k,:) = strfield(3:5);
        else
            coord(k,:) = strfield(4:6);
        end
        k = k + 1;
    end
    
end

fclose(fid);

status = 0;
if( length(coord) > 0 & ~isConverged )
    fprintf('Warning: Optimization is not finished yet in file ''%s''!!!\n', inputFileName);
    status = 1;
end
