%% read images
pic_1 = imread('room1.jpg');
pic_1 = rgb2gray(pic_1);
imshow(pic_1);

pic_2 = imread('room2.jpg');
pic_2 = rgb2gray(pic_2);


%% SURF
pts_1 = detectSURFFeatures(pic_1);
pts_2 = detectSURFFeatures(pic_2);

[features_1,  valid_pts_1]  = extractFeatures(pic_1,  pts_1);
[features_2,  valid_pts_2]  = extractFeatures(pic_2,  pts_2);

indexPairs = matchFeatures(features_1, features_2);

matched_1  = valid_pts_1(indexPairs(:,1));
matched_2 = valid_pts_2(indexPairs(:,2));

M_1 = getNorm(matched_1);
M_2 = getNorm(matched_2);

% normalize
matched_1_applied = applyM(M_1, matched_1.Location);
matched_2_applied = applyM(M_2, matched_2.Location);

%% RANSAC
h_final = [];
c_final = 0;
% maximum iteration: 100
for i=1:100
    idx = randperm(length(matched_1), 4);
    matched_1_set = matched_1.Location(idx, :);
    matched_2_set = matched_2.Location(idx, :);
    
    matched_1_set_applied = matched_1_applied(idx, :);
    matched_2_set_applied = matched_2_applied(idx, :);
    
    mat_a = generateA(4, matched_1_set_applied, matched_2_set_applied);
    h = getHNormalized(mat_a, M_1, M_2);
    
    tmp = applyM(h, matched_1.Location);
    dist_val = (tmp - matched_2.Location) .^ 2;
    dist_val = sqrt(sum(dist_val, 2));
    dist = dist_val < (sqrt(size(pic_1,1)^2 + size(pic_2,2)^2) / 100);
    c_sum = sum(dist);
    
    if c_sum > 10
        mat_a = generateA(c_sum, matched_1_applied(dist, :), matched_2_applied(dist, :));
        h = getHNormalized(mat_a, M_1, M_2);
        
        if c_sum > c_final
            c_final = c_sum;
            h_final = h;
        end
    end
end

%% Image Stitching

x_offset = round(size(pic_1,1) * 0.8);
y_offset = round(size(pic_1,2) * 0.8);

fill = zeros(size(pic_1));
pic_1_red = cat(3, pic_1, fill, fill);
im_out = padarray(pic_1_red, [y_offset x_offset], 0, 'both');

for y=1:size(im_out, 1)
    for x=1:size(im_out, 2)
       pt_1 = [x-x_offset y-y_offset];
       pt_2 = round(applyM(h_final, pt_1))';
       
       if pt_2(1) > 0 && pt_2(2) > 0
           if ~(pt_2(1) > size(pic_2, 2)) && ~(pt_2(2) > size(pic_2, 1))
               % inside second image
               im_out(y, x, 2) = pic_2(pt_2(2), pt_2(1));
               im_out(y, x, 3) = pic_2(pt_2(2), pt_2(1));
           end
       end
    end
end

figure;
imshow(im_out);


%% Helper Functions
function norm = getNorm(match)
    mat = match.Location;
    n = size(mat, 1);
    mat_mean = mean(mat);
    new_mat = (mat - mat_mean) .^ 2;
    mat_sum = sum(new_mat, 'all') / (2*n);
    std = sqrt(mat_sum);
    T = [1 0 mat_mean(1); 0 1 mat_mean(2); 0 0 1];
    S = [1/std 0 0; 0 1/std 0; 0 0 1];
    norm = S * T;
end

function A = generateA(size, matched_1, matched_2)
    % Directly from page 23, lecture 17
    A = zeros(2 * size, 9);
    for i=1:size
        x_pos = matched_1(i, 1);
        y_pos = matched_1(i, 2);
        x_pos_2 = matched_2(i, 1);
        y_pos_2 = matched_2(i, 2);
        A(2*i - 1, :) = [x_pos y_pos 1 0 0 0 -x_pos_2*x_pos -x_pos_2*y_pos -x_pos_2];
        A(2*i, :) = [0 0 0 x_pos y_pos 1 -y_pos_2*x_pos -y_pos_2*y_pos -y_pos_2]; 
    end
end

function H = getHNormalized(mat_in, mat_1, mat_2)
    [U, S, V] = svd(mat_in);
    sol = V(:,end);
    H = mat_2 \ (reshape(sol,[],3)') * mat_1;
end

function out = applyM(M, in_pts)
    in_pts = [in_pts ones(size(in_pts, 1), 1)]';
    in_pts = M * in_pts;
    out = [in_pts(1,:) ./ in_pts(3,:); in_pts(2,:) ./ in_pts(3,:)]';  
end
