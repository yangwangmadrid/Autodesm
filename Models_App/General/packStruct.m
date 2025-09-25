%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

function S = packStruct(varargin)
    S = struct();
    for k = 1:nargin
        name = inputname(k);
        if isempty(name)
            error("All inputs must be named variables.");
        end
        S.(name) = varargin{k};
    end
end