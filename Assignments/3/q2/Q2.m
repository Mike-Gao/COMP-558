%% Loading Images

im_1 = imread('bookshelf1.jpg');
im_1 = rgb2gray(im_1);

im_2 = imread('bookshelf2.jpg');
im_2 = rgb2gray(im_2);

figure;
subplot(1,2,1);
imshow(im_1);
title("Image 1");

subplot(1,2,2); 
imshow(im_2);
title("Image 2");

%% SURF

pts_1 = detectSURFFeatures(im_1);
pts_2 = detectSURFFeatures(im_2);

[features_1,  valid_pts_1]  = extractFeatures(im_1,  pts_1);
[features_2,  valid_pts_2]  = extractFeatures(im_2,  pts_2);

indexPairs = matchFeatures(features_1, features_2);

matched_1  = valid_pts_1(indexPairs(:,1));
matched_2 = valid_pts_2(indexPairs(:,2));


%% Estimate Fundamental Matrix

[F,inliersIndex] = estimateFundamentalMatrix(matched_1,matched_2, 'Method', 'RANSAC');

inliers_1 = matched_1(inliersIndex, :);
inliers_2 = matched_2(inliersIndex, :);

figure;
showMatchedFeatures(im_1, im_2, inliers_1, inliers_2, 'montage');
title('Inlier Point Matches');

%% Estimate uncalibrated rectification

[T1,T2] = estimateUncalibratedRectification(F,inliers_1,inliers_2,size(im_2));

[im_1_rectified, im_2_rectified] = rectifyStereoImages(im_1, im_2, T1, T2);

figure;
imshow(stereoAnaglyph(im_1_rectified, im_2_rectified));
title('Rectified Image');

%% Disparity

disparityMap = disparitySGM(im_1_rectified, im_2_rectified, 'DisparityRange', [-56 56], 'UniquenessThreshold', 6);

figure;
imshow(disparityMap, [-56 56]);
title('Disparity Map');
colormap jet;
colorbar;

%% Epipole

num_lines = size(inliers_1, 1);
L1 = zeros(num_lines, 3);
L2 = zeros(num_lines, 3);
epi_line_start_1 = zeros(num_lines, 3); 
epi_line_start_2 = zeros(num_lines, 3);
epi_line_end_1 = zeros(num_lines, 3); 
epi_line_end_2 = zeros(num_lines, 3);
in_1 = zeros(num_lines, 2);
in_2 = zeros(num_lines, 2);
for i = 1:num_lines
    in_1(i,:)= inliers_1.Location(i,:);
    in_2(i,:)= inliers_2.Location(i,:);
    
    L1(i,:) = [inliers_2.Location(i,:), 1] * F;
    L2(i,:) = F * [inliers_1.Location(i,:), 1]';
    epi_line_start_1(i, :) = [0, -L1(i,3)/L1(i,2), 1];
    epi_line_end_1(i, :) = [size(im_1,2), -size(im_1,2)*L1(i,1)-L1(i,3)/L1(i,2), 1];
    epi_line_start_2(i, :) = [0, -L2(i,3)/L2(i,2), 1];
    epi_line_end_2(i, :) = [size(im_2,2), -size(im_2,2)*L2(i,1)-L2(i,3)/L2(i,2), 1];
end

[~, ~, V1] = svd(L1);
[~, ~, V2] = svd(L2);

epipole_1 = V1(:,end);
epipole_1_divided_by_w = [epipole_1(1)/epipole_1(3), epipole_1(2)/epipole_1(3)];

epipole_2 = V2(:,end);
epipole_2_divided_by_w = [epipole_2(1)/epipole_2(3), epipole_1(2)/epipole_2(3)];


figure;
imshow(im_1);
title("Epipole in Image 1 is x= " + epipole_2_divided_by_w(1) + " y = " + epipole_2_divided_by_w(2));

% limit to the line print out to 10
line_lim = 10;
color = hsv(line_lim);

for i=1:line_lim
    viscircles(in_1(i,:), 10, 'Color', color(i,:));
    l = line([epi_line_start_1(i,1), epi_line_end_1(i,1)], [epi_line_start_1(i,2), epi_line_end_1(i,2)]);
    l.Color = color(i,:);
end


figure;
imshow(im_2);
title("Epipole in Image 2 is x= " + epipole_2_divided_by_w(1) + " y = " + epipole_2_divided_by_w(2));

for i=1:line_lim
    viscircles(in_2(i,:), 10, 'Color', color(i,:));
    l = line([epi_line_start_2(i,1), epi_line_end_2(i,1)], [epi_line_start_2(i,2), epi_line_end_2(i,2)]);
    l.Color = color(i,:);
end




