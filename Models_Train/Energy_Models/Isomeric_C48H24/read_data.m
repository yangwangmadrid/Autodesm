%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

%

function data = read_data()

dir_data = '../../../Datasets/Benzenoids_C48H24';
inp = sprintf( '%s/DATA_complete.dat', dir_data );
data = readData( inp );
fprintf( 'A total of %i PAH molecules\n', length(data.name_arr) )

end

function data = readData( inp )

name_arr = {};
plan_arr = [];
EDFT_arr = [];
EHMO_arr = [];
NC_arr = [];
NH_arr = [];
distHH_arr = {};
CCBondLen_arr = {};
CHBondLen_arr = {};

fid = fopen( inp );
%k = 0;
while ~feof(fid)
    tline = fgetl(fid);
    s = strsplit( strtrim(tline) );
    % Name:
    name_arr{end+1,1} = s{1};
    % Planarity:
    plan_arr(end+1,1) = str2double( s{2} );
    % E_DFT:
    EDFT_arr(end+1,1) = str2double( s{3} );
    % E_HMO:
    EHMO_arr(end+1,1) = str2double( s{4} );
    % Number of C atoms:
    NC_arr(end+1,1) = str2double( s{5} );
    % Number of H atoms:
    NH_arr(end+1,1) = str2double( s{6} );
    % distHH:
    nHH = str2double( s{7} );
    distHH = zeros(nHH,1);
    for j = 1 : nHH
        distHH(j) = str2double( s{7+j} );
    end
    distHH_arr{end+1,1} = distHH;
    % CCBondLen:
    nCCBond = str2double( s{nHH+8} );
    CCBondLen = zeros(nCCBond,1);
    for j = 1 : nCCBond
        CCBondLen(j) = str2double( s{nHH+8+j} );
    end
    CCBondLen_arr{end+1,1} = CCBondLen;
    % CHBondLen:
    nCHBond = str2double( s{nHH+nCCBond+9} );
    assert( nCHBond == NH_arr(end) )
    CHBondLen = zeros(nCHBond,1);
    for j = 1 : nCHBond
        CHBondLen(j) = str2double( s{nHH+nCCBond+9+j} );
    end
    CHBondLen_arr{end+1,1} = CHBondLen;


    % k = k + 1;
    % if k > 5
    %     break
    % end
end
fclose( fid );

data.name_arr = name_arr;
data.plan_arr = plan_arr;
data.EDFT_arr = EDFT_arr;
data.EHMO_arr = EHMO_arr;
data.NC_arr = NC_arr;
data.NH_arr = NH_arr;
data.distHH_arr = distHH_arr;
data.CCBondLen_arr = CCBondLen_arr;
data.CHBondLen_arr = CHBondLen_arr;

end


function data = mergeData( data_arr )

data.name_arr = {};
data.plan_arr = [];
data.EDFT_arr = [];
data.EHMO_arr = [];
data.NC_arr = [];
data.NH_arr = [];
data.distHH_arr = {};
data.CCBondLen_arr = {};
data.CHBondLen_arr = {};

fprintf( 'Merging data ...\n' )
for j = 1 : length( data_arr )
    N = length( data_arr(j).name_arr );
    for k = 1 : N
        data.name_arr{end+1,1} = data_arr(j).name_arr{k};
        data.plan_arr(end+1,1) = data_arr(j).plan_arr(k);
        data.EDFT_arr(end+1,1) = data_arr(j).EDFT_arr(k);
        data.EHMO_arr(end+1,1) = data_arr(j).EHMO_arr(k);
        data.NC_arr(end+1,1) = data_arr(j).NC_arr(k);
        data.NH_arr(end+1,1) = data_arr(j).NH_arr(k);
        data.distHH_arr{end+1,1} = data_arr(j).distHH_arr{k};
        data.CCBondLen_arr{end+1,1} = data_arr(j).CCBondLen_arr{k};
        data.CHBondLen_arr{end+1,1} = data_arr(j).CHBondLen_arr{k};
    end
end

end