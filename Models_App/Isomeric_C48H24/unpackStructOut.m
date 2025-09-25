%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

function varargout = unpackStructOut(S, varargin)
    % Return struct fields as ordered outputs
    if nargin < 2
        fields = fieldnames(S);
    else
        fields = varargin;
    end

    n = numel(fields);
    varargout = cell(1, n);
    for k = 1:n
        varargout{k} = S.(fields{k});
    end
end
