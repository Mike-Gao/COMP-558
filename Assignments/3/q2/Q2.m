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

disparityMap = disparitySGM(im_1_rectified, im_2_rectified, 'DisparityRange', [-55 55], 'UniquenessThreshold', 6);

figure;
imshow(disparityMap, [-55 55]);
title('Disparity Map');
colormap jet;
colorbar;

%% Epipole




