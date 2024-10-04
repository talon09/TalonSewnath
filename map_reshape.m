function [reshaped_map] = map_reshape(xyzI)
    % INPUT: n by 4 matrix of x y z and intensity
    % OUTPUT: 7 by 7 by 4 
    reshaped_map = zeros(7,7,4);
    sorted = sortrows(xyzI, 4, "descend"); 
    batched = sorted(1:49,:);
    batched_sorted = sortrows(batched, [1,3,2,4]); % Sort on width, depth, height and intensity

    for feature_idx = 1:4
        feature_column = batched_sorted(:, feature_idx);
        % Reshape the feature column into an 7x7 grid
        grid = reshape(feature_column, [7, 7]);
        reshaped_map(:, :, feature_idx) = grid;
    end

end