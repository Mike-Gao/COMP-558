% Load Image %
im_syn = rgb2gray(imread('synthetic.jpg'));
im_nat = rgb2gray(imread('mtl.jpg'));
% Q7 - resized image
% im_syn = rgb2gray(imresize(imread('synthetic.jpg'),0.8));
gss_syn = getgss(im_syn);
gss_nat = getgss(im_nat);

dog_syn = getdog(gss_syn);
dog_nat = getdog(gss_nat);

% Q1: Gaussian Scale Space
im = rgb2gray(imread('synthetic.jpg'));
figure
title("synthetic");
for i = 1:16
    sig = getsig(i+1);
    gauss = fspecial('gaussian',floor(sig*2),sig);
    newim = imfilter(im,gauss);
    subplot(4,4,i)
    imshow(newim)
    title(sprintf("Gaussian %i: Sigma value %f",i+1,sig));
end 

nat = rgb2gray(imread('mtl.jpg'));
figure
title("natural");
for i = 1:16
    sig = getsig(i+1);
    gauss = fspecial('gaussian',floor(sig*2),sig);
    newnat = imfilter(nat,gauss);
    subplot(4,4,i)
    imshow(newnat)
    title(sprintf("Gaussian %i: Sigma value %f",i+1,sig));
end 


% Q2: Harris-Stevens 
figure
for i = 1:16
    sig = getsig(i+1);
    subplot(4,4,i)
    plotimage(harris(2*sig, squeeze(gss_syn(:,:,(i+1)))), 'jet', sprintf('Sigma: %f', sig))
end

figure
for i = 1:16
    sig = getsig(i+1);
    subplot(4,4,i)
    plotimage(harris(2*sig, squeeze(gss_nat(:,:,(i+1)))), 'jet', sprintf('Sigma: %f', sig))
end



% Q3: Difference of Gaussians                   
figure
for i = 1:16
    subplot(4,4,i)
    plotimage(squeeze(dog_syn(:,:,i)), 'jet', sprintf('Layer: %d', i))
end

figure
for i = 1:16
    subplot(4,4,i)
    plotimage(squeeze(dog_nat(:,:,i)), 'jet', sprintf('Layer: %d', i))
end


% Q4: SIFT Keypoint detection
figure
imagesc(im_syn);
colormap(gray);
keypts_syn = getkeypts(dog_syn);
for i=1:size(keypts_syn,1)
     viscircles([keypts_syn(i,1), keypts_syn(i,2)], getsig(keypts_syn(i,3)));
end

% figure
% imagesc(im_nat);
% colormap(gray);
% keypts_nat = getkeypts(dog_nat);
% for i=1:size(keypts_nat,1)
%      viscircles([keypts_nat(i,1), keypts_nat(i,2)], keypts_nat(i,3));
% end

% Q5: Hessian constraint
figure
imagesc(im_syn);
colormap(gray);
xfilter = [1 0 -1; 1 0 -1; 1 0 -1]/6;
yfilter = [1 0 -1; 1 0 -1; 1 0 -1]'/6;
r = 10;
threshold = (r+1)^2/r;
keypoints = keypts_syn;
keypoints_new = zeros(length(keypoints),4);
cur = 1;
for i = 1:size(keypoints,1)
    x = keypoints(i,1);
    y = keypoints(i,2);
    l = keypoints(i,3);
    neighbors = dog_syn(y-2:y+2, x-2:x+2, l);
    dx = conv2(neighbors,xfilter,'valid');
    dy = conv2(neighbors,yfilter,'valid');
    dxx = conv2(dx, xfilter, 'valid');
    dyy = conv2(dy, yfilter, 'valid');
    dxy = conv2(dy, yfilter, 'valid');
    hessian = [dxx, dxy; dxy, dyy];
    if (trace(hessian)^2/det(hessian) < threshold && det(hessian) > 0)
        keypoints_new(cur,:) = keypoints(i,:);
        cur = cur + 1;
        viscircles([x, y], getsig(l), 'Color','r');
    else
        viscircles([x, y], getsig(l), 'Color','b');
    end
end
keypts_syn = keypoints_new(1:cur-1,:);


% figure
% imagesc(im_nat);
% colormap(gray);
% xfilter = [1 0 -1; 1 0 -1; 1 0 -1]/6;
% yfilter = [1 0 -1; 1 0 -1; 1 0 -1]'/6;
% r = 10;
% threshold = (r+1)^2/r;
% keypoints = keypts_nat;
% keypoints_new = zeros(length(keypoints),4);
% cur = 1;
% for i = 1:size(keypoints,1)
%     x = keypoints(i,1);
%     y = keypoints(i,2);
%     l = keypoints(i,3);
%     neighbors = dog_nat(y-2:y+2, x-2:x+2, l);
%     dx = conv2(neighbors,xfilter,'valid');
%     dy = conv2(neighbors,yfilter,'valid');
%     dxx = conv2(dx, xfilter, 'valid');
%     dyy = conv2(dy, yfilter, 'valid');
%     dxy = conv2(dy, yfilter, 'valid');
%     hessian = [dxx, dxy; dxy, dyy];
%     if (trace(hessian)^2/det(hessian) < threshold && det(hessian) > 0)
%         keypoints_new(cur,:) = keypoints(i,:);
%         cur = cur + 1;
%         viscircles([x, y], getsig(l), 'Color','r');
%     else
%         viscircles([x, y], getsig(l), 'Color','b');
%     end
% end

% Q6: SIFT feature dominant orientation
figure;
imagesc(im_syn);
colormap(gray);
binhists = zeros(300, 36);
sift_feature_dominant = zeros(1000,4);
keypoints = keypts_syn;
for i = 1:size(keypoints,1)
    x = keypoints(i,1);
    y = keypoints(i,2);
    l = keypoints(i,3);
    sig_round = round(getsig(l));
    neighbors = gss_syn(x-sig_round-1:x+sig_round+1, y-sig_round-1:y+sig_round+1,l);
    dx = conv2(neighbors,xfilter,'valid');
    dy = conv2(neighbors,yfilter,'valid');
    mag = sqrt(dx.*dx + dy.*dy);
    gauss = fspecial('gaussian', floor(2*getsig(l)), 1.5*getsig(l));
    mag_weighted = mag .* gauss;
    dir = atan2d(-dx, dy);
    dir(dir<0) = dir(dir<0) + 360;
    bin = round(dir/10) + 1;
    bin(bin==37) = 1;
    bin_hist = zeros(1,36);
    for j = 1:size(mag_weighted,1)
        bin_hist(bin(j)) = bin_hist(bin(j)) + mag_weighted(j);
    end
    binhists(pts,:) = bin_hist;
    pos = pos + 1;
    large_vals = binhist( binhist > 0.8 * max(binhist) )';
    for k=1:length(large_vals)
        idx = find(binhist == larges(k));
        ang = (idx-1)*10;
        sift_feature_dominant(cur,:) = [x y l ang];
        cur = cur + 1;
        delta_x = 2*getsig(l)*cos(-deg2rad(ang));
        delta_y = 2*getsig(l)*sin(-deg2rad(ang));
        line([y y+delta_y], [x x+delta_x],'Color', 'r', 'LineWidth',2);
    end
end

    
% Q7 - Scale Invariance
% Change im_syn to rgb2gray(imresize(imread('synthetic.jpg'),0.8)), see
% above.



% Helper functions                     
function sig = getsig(i)
    k = floor(i/4);
    m = mod(i,4);
    sig = 2^(k+m/4);
end

function gss = getgss(im)
   gss = zeros(size(im,1), size(im,2),17, 'uint8');
   for i = 1:17
       sig = getsig(i);
       gauss = fspecial('gaussian', floor(sig*2), sig);
       new = imfilter(im, gauss);
       gss(:,:,i) = new;
   end
end

function hs = harris(sig,im)
    [Gx,Gy] = imgradientxy(im);
    GxSq = Gx.*Gx;
    GySq = Gy.*Gy;
    GxGy = Gx.*Gy;

    outerGauss = fspecial('gaussian',floor(sig),sig);

    MxSq = conv2(GxSq,outerGauss);
    MySq = conv2(GySq,outerGauss);
    MxMy = conv2(GxGy,outerGauss);

    t = size(im,1);
    s = size(im,2);
    hs = zeros(t,s);
    for x = 1:s
        for y = 1:t
            M = [MxSq(y,x) MxMy(y,x) ; MxMy(y,x) MySq(y,x)];
            hs(y,x) = det(M) - 0.1 * trace(M)^2;
        end
    end
end

function dog = getdog(gss)
dog = zeros(size(gss,1), size(gss,2), size(gss,2), 'uint8');
    for i = 1:16
        dog(:,:,i) = gss(:,:,i) - gss(:,:,i+1);
    end
end

function pts = getkeypts(dog)
h = size(dog, 1);
w = size(dog, 2);
lim = 20;
pts = zeros(10000,4);
cur = 1;
for i = 2:15
    sig = getsig(i);
    offset = ceil(2*sig)+1;
    for x = offset: w - offset
        for y = offset: h - offset
               if abs(dog(y,x,i)) > lim
                   neighbors = dog(y-1:y:y+1,x-1:x:x+1,i-1:i+1);
                   max_neighbor = max(neighbors, [], 'all');
                   min_neighbor = min(neighbors, [], 'all');
                   if (dog(y,x,i) == max_neighbor || dog(y,x,i) == min_neighbor)
                       pts(cur,1:3) = [x y i];
                       cur = cur + 1;
                   end
               end
        end
    end
end
pts = pts(pts(:,3) ~= 0,:); 
end

function plotimage(im,map,t)
    imagesc(im)
    colormap(map)
    colorbar
    title(t)
end
    
