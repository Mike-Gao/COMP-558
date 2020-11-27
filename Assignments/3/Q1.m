
pic_1 = imread('yard1.jpg');
pic_1 = rgb2gray(pic_1);
imshow(pic_1);

pic_2 = imread('yard2.jpg');
pic_2 = rgb2gray(pic_2);
imshow(pic_2);

pts_1 = detectSURFFeatures(pic_1);
pts_2 = detectSURFFeatures(pic_2);

[features_1,  valid_pts_1]  = extractFeatures(pic_1,  pts_1);
[features_2,  valid_pts_2]  = extractFeatures(pic_2,  pts_2);

indexPairs = matchFeatures(features_1, features_2);

matched_1  = valid_pts_1(indexPairs(:,1));
matched_2 = valid_pts_2(indexPairs(:,2));

figure;
showMatchedFeatures(pic_1,pic_2,matched_1,matched_2);


idx_random_pair = randperm(length(indexPairs), 4);

out = getNorm(matched_1);

function norm = getNorm(match)
    mat = match.Location;
    n = size(mat, 1);
    mat_mean = mean(mat);
    new_mat = (mat - mat_mean) .^ 2;
    mat_sum = sum(new_mat) / (2*n);
    std = sqrt(mat_sum);
    T = [1 0 mat_mean(1); 0 1 mat_mean(2); 0 0 1];
    S = [1/std 0 0; 0 1/std 0; 0 0 1];
    norm = S * T;
end

function A = generateA(size, matched_1, matched_2)
    % Directly from page 23, lecture 17
    A = zeros(2 * size, 9);
    for i=1:size
        x_pos = matched_1.Location(i, 1);
        y_pos = matched_1.Location(i, 2);
        x_pos_2 = matched_2.Location(i, 1);
        y_pos_2 = matched_2.Location(i, 2);
        A(2*i - 1, :) = [x_pos y_pos 1 0 0 0 -x_pos_2*x_pos -x_pos_2*y_pos -x_pos_2];
        A(2*i, :) = [0 0 0 x_pos y_pos 1 -y_pos_2*x_pos -y_pos_2*y_pos -y_pos_2]; 
    end
end