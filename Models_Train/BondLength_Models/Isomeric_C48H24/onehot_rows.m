%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

function [onehot, labels] = onehot_rows(A)
    % Ensure A is a column vector
    A = A(:);

    % Get sorted unique values
    labels = unique(A);
    nLabels = numel(labels);
    nInstances = numel(A);

    % Map each A(i) to its index in labels
    [~, idx] = ismember(A, labels);

    % Create one-hot matrix (rows = instances, columns = labels)
    onehot = zeros(nInstances, nLabels);
    for i = 1:nInstances
        onehot(i, idx(i)) = 1;
    end
end
