%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

% UNPACKSTRUCT Unpack structure fields to variables:
function varargout = unpackStruct(s)

fields = fieldnames(s);
varargout = cell(1, nargout);
for i = 1:nargout
    varargout{i} = s.(fields{i});
end

end