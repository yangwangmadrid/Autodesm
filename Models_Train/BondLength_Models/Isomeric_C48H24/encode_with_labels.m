%% Source codes: Autodesmotic reactions for strain energy evaluation
%% Author: Yang Wang (yangwang@yzu.edu.cn)
%% Last modified: September 25, 2025
%% License: For academic and non-commercial use only

function onehot = encode_with_labels(labels, a)
    % Ensure a is a column vector
    a = a(:);
    nInstances = numel(a);
    nLabels = numel(labels);

    % Map each value in 'a' to index in 'labels'
    [~, idx] = ismember(a, labels);

    % Warn if any value in 'a' is not in 'labels'
    if any(idx == 0)
        setdiff( a, intersect( labels, a ) )
        error('Some values in ''a'' are not in ''labels'' and will be encoded as all zeros.');
    end

    % Initialize one-hot matrix
    onehot = zeros(nInstances, nLabels);

    % Fill one-hot matrix where label match is found
    for i = 1:nInstances
        if idx(i) > 0
            onehot(i, idx(i)) = 1;
        end
    end
end