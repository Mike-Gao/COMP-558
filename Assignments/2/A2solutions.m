%%
% Assignment 2 solutions

clear all

chooseImage

display('for the resizing question -- uncomment line 9 in the code')
I = uint8(imresize(I, .6));

Nx = size(I,2);
Ny = size(I,1);

figure
imagesc(I)
xlabel('x','fontsize',14);
ylabel('y','fontsize',14);
colormap gray

%%
% Q1:   compute a Gaussian scale space

close all
%  Scales used
sigmas = power(2, 1 : 0.25 : 5);

numSigmas = size(sigmas,2);
gaussianScaleSpace = zeros(numSigmas, Ny, Nx );

%   Compute the Gaussian scale space.

for sigmaCt = 1:numSigmas
    sig = sigmas(sigmaCt);
    %  blur the image with a Gaussian of standard deviation sig, and make
    %  the radius of this Gaussian be 3*sig
    g = fspecial('gaussian', 2*ceil(3*sig)+1, sig);
    gaussianScaleSpace(sigmaCt,:,:) = conv2(I,g, 'same');
end
figure
colormap gray
subplotScalespace( gaussianScaleSpace(2:end,:,:), 'Gaussian scale space');

%%
% Compute the Harris corner scale space

HarrisScaleSpace = zeros(numSigmas, Ny, Nx );
IxScaleSpace = zeros(numSigmas, Ny, Nx );
IyScaleSpace = zeros(numSigmas, Ny, Nx );

for sigmaCt = 1:numSigmas
    %  compute gradient at each scale
    Ix = conv2(squeeze(gaussianScaleSpace(sigmaCt,:,:)), [1, 0, -1], 'same');
    Iy = conv2(squeeze(gaussianScaleSpace(sigmaCt,:,:)), [1, 0, -1]', 'same');
    
    IxScaleSpace(sigmaCt,:,:) = Ix;
    IyScaleSpace(sigmaCt,:,:) = Iy;
    
    %  compute the second moment matrix by blurring the structure tensor
    %  at each scale with a Gaussian that is twice as large as the inner scale Gaussian
    
    M_11 = Ix .* Ix;
    M_22 = Iy .* Iy;
    M_cross = Ix .* Iy;
    
    sig = sigmas(sigmaCt);
    %  The window for the 2nd moment matrix should have radius at least twice the
    %  scale (2*sig) at this level.
    %  Note this will not give meaningful results near the image boundary.
    
    g = fspecial('gaussian', 2*ceil(2*sig)+1, 2*sig);
    M_11 = conv2(M_11, g,'same');
    M_22 = conv2(M_22, g,'same');
    M_cross = conv2(M_cross, g, 'same');
    
    detM = M_11 .* M_22 - M_cross .* M_cross;
    trM  = M_11 + M_22;
    HarrisScaleSpace(sigmaCt,:,:) = detM - 0.1 * trM.*trM;
end

figure
colormap jet
subplotScalespace( HarrisScaleSpace(2:end,:,:), 'Harris-Stevens corner scale space');
%%
% Compute the DOG scale space  
%  16 layers -- only find keypoints in layers 2 to 15
%  Note I wrote the code to be a bit more general than that, allowing for
%  5x5x5 neighborhoods too.

DoGI = zeros(numSigmas-1,Ny, Nx );
for sigmaCt = 2:numSigmas
    sigmaD = sigmas(sigmaCt);
    DoGI(sigmaCt,:,:) = gaussianScaleSpace(sigmaCt,:,:) - gaussianScaleSpace(sigmaCt-1,:,:);
end
figure
colormap jet
subplotScalespace( DoGI(2:end,:,:), 'Difference of Gaussians scale space');

%%   Find the keypoints

%Neighborhood is scale space is a cube of width  2*radius+1, 3
radius = 1;   % 2 ->  5x5x5 cube
% 1 ->  3x3x3 cube

% make a big array to hold the keypoints
%  I only do this because Matlab doesn't like variable size lists.
%  Once I've computed all the keypoints, I'll resize the list.
maxNumKeypoints = 10000;

%  keypoint is  (sigma, x, y, theta)
keypoints = zeros(maxNumKeypoints,4);
keypointCt = 0;

%  Keep track of keypoint that fail the Hessian test so that I can color
%  them blue

HessFailedKeypoints = zeros(maxNumKeypoints,4);
HessFailedCt = 0;

%  The keypoints need to be minima or maxima, but we want to have DOG
%  values sufficiently different from 0.  The choice of threshold here is
%  arbitrary.

threshAbsDoG =  1.0*std(DoGI(:));

%  We are looking for local maxima and minima in a 3x3x3 neighborhood (radius = 1),
%  so we cannot use the first or last scale.

for sigmaCt = 1+radius:numSigmas-radius-1
    
    %   Only consider keypoints whose position (x,y) is sufficiently far away
    %   from the image boundary.
    %   The orientation selection step will check gradient vectors in a
    %   larger neighborhood of keypoint, so set the boundary width to be
    %   sufficiently large:  3 sigma.
    
    boundaryWidth = max(radius+1, ceil( sigmas(sigmaCt)*3 ));
    
    for i = boundaryWidth+1:Nx-boundaryWidth-1
        for j = boundaryWidth+1:Ny-boundaryWidth-1
            %  only consider keypoints for which the DoGI is above some
            %  threshold
            if (abs(DoGI(sigmaCt,j,i)) > threshAbsDoG)
                neigh = reshape( DoGI(sigmaCt-radius:sigmaCt+radius,j-radius:j+radius,i-radius:i+radius), [1 power(2*radius+1,3) ]);
                neigh  = [neigh(1:floor(power(2*radius+1,3)/2)) neigh(ceil(power(2*radius+1,3)/2 + 1): power(2*radius+1,3))];
                
                val = DoGI(sigmaCt,j,i);
                
                if (val > max(neigh)) || (val < min(neigh))
                    
                    %  Lowe's suggests a condition on the Hessian, namely
                    %    ensure tr(Hessian)^2 / det(Hessian) < (r+1)^2/r                  
                    Hess_xx = DoGI(sigmaCt, j, i+1) - 2*DoGI(sigmaCt, j, i) + DoGI(sigmaCt, j, i-1);
                    Hess_yy = DoGI(sigmaCt, j+1, i) - 2*DoGI(sigmaCt, j, i) + DoGI(sigmaCt, j-1, i);
                    Hess_xy = ((DoGI(sigmaCt, j+1, i+1) -  DoGI(sigmaCt, j-1, i+1))/2 - ...
                        (DoGI(sigmaCt, j+1, i-1) -  DoGI(sigmaCt, j-1, i-1))/2) / 2;
                    % we take r=10, as per the paper, thus (r+1)^2/r=12.1
                    testHess = (power(Hess_xx + Hess_yy,2) / (Hess_xx * Hess_yy - Hess_xy * Hess_xy) < 12.1) ...
                        &&  (Hess_xx * Hess_yy - Hess_xy * Hess_xy  >= 0);
                    
                    %   ASIDE:  it can happen that the determinant is less than 0, even though the point is local maximum or minimum.
                    %   The Hessian doesn't care about a linear  fit.   So when you subtract off the linear fit,
                    %   the point might not be a local maximum or minimum anymore, and so determinant can be negative.
                    %   This seems to give lousy keypoints.
                   
                    if testHess
                        %  The keypoint candidate passes the Hessian condition
                        keypointCt = keypointCt+1;
                        %  We will compute the dominant orientation later.
                        keypoints(keypointCt,:) = [sigmaCt, i, j, 0];
                    else
                        %  the keypoint candidate fails the Hessian condition
                        HessFailedCt = HessFailedCt + 1;
                        HessFailedKeypoints(HessFailedCt, :) = [sigmaCt, i, j, 0];
                    end
                end
            end
        end
    end
end

%%  Compute the dominant orientation at each keypoint.

% For this, we will look at the gradient vectors in a Gaussian neighborhood
% of radius of 3*sigma
%
% For each keypoint, we compute the orientation histogram and find where the maximum value occurs.
% If there is another orientation at which the histogram is within some
% tolerance threshold of the max,  then we'll include this keypoint also
% (by appending to the end of the keypoint list).

extraKeypointCt = keypointCt;
peak_hist_thresh = 0.9 ; 

for ct = 1:keypointCt
    sigIdx = keypoints(ct,1);
    x0 = keypoints(ct,2);
    y0 = keypoints(ct,3);
    sig0  = sigmas( sigIdx );
    sigWindow = sig0*1.5;     
    %Lowe (2004) says to use a Gaussian window that is 1.5 times scale of
    %keypoint.   I do so, but only take up to 2 standard deviations.
    
    lowY = floor(y0 - 2*sigWindow);
    highY = ceil(y0 + 2*sigWindow);
    lowX = floor(x0 - 2*sigWindow);
    highX = ceil(x0 + 2*sigWindow);

    %  for plot comparison only
    localI = I( lowY:highY, lowX:highX);
    localIGaussScaleSpace = squeeze(gaussianScaleSpace( sigIdx , lowY:highY, lowX:highX));
    Ix = squeeze(IxScaleSpace(sigIdx,lowY:highY, lowX:highX));
    Iy = squeeze(IyScaleSpace(sigIdx,lowY:highY, lowX:highX));
    width = size(Ix,1);
    g = fspecial('gaussian',width, sigWindow);
    theta = atan2(Iy, Ix);   %  theta goes from -pi to pi i.e. -180 deg to 180 deg
    
    % gaussian weight the gradients
    Ix = Ix .* g;
    Iy = Iy .* g;
    
    %Make the first bin have theta =0 at its center
    
    % theta = 0 will map to bin 18
    %thetaIndex = ceil(  mod((theta/pi + 1)*18 - .5, 36));
    
    % theta = 0 will map (-5,5] deg to bin 17.
    % Later I will circshift the histogram so theta=0 will map to bin 1
    thetaIndex = ceil(  mod((theta/pi + 1)*18 + .5, 36));
    
    mag =  sqrt(Ix .* Ix + Iy .* Iy) ;
    orientationHist = zeros(36,1);
    for j = 1:width
        for i= 1:width
            orientationHist(thetaIndex(j,i)) = orientationHist(thetaIndex(j,i)) + mag(j,i) ;
        end
    end
    %  This will map theta=0 to bin 1 out of 36.   
    thetaIndex = circshift(thetaIndex, 18); 
    
    % find the maximum value
    maxHistValue = max(orientationHist(:));
    idx = find(maxHistValue == orientationHist);
    if length(idx) > 1
        idx = idx(1);
    end
    % idx is in [1, 36]
    % This is the theta index of the keypoint.
    keypoints(ct,4) =  (idx - 1)/18.0*pi;   %  gives a theta value in [0, 2 pi)
    
    % Here I want to add in code that will give multiple keypoints for 
    % orientation peaks greater than peak_thres of max. 
    
    for k = 1:36
        if (k ~= idx) && (orientationHist(k) >= peak_hist_thresh * maxHistValue)
            extraKeypointCt = extraKeypointCt + 1;
            for m = 1:3
                keypoints(extraKeypointCt,m) = keypoints(ct,m);
            end
            keypoints(extraKeypointCt,4) = (k - 1)/18.0*pi; %  gives a theta value in [0, 2 pi)
        end
    end
    
    %  2x2 plot of local image patch and corresponding patch in scale space
    %  along with vector field and orientation histogram
    %
    %  I used this for the lecture slides.   But it is not part of the
    %  assignment.
    
    if (0)
        figure
        subplot(2,2,1)
        imagesc( localI  );  axis square; colormap gray
        title('original image (cropped)')
        subplot(2,2,2)
        imagesc( localIGaussScaleSpace  );  axis square; colormap gray
        title('image at scale of keypoint')
        subplot(2,2,3)
        %  warn students that this is nasty
        %[x,y] = meshgrid(1:size(localI,1),1:size(localI,1));
        [x,y] = meshgrid(2:size(localI,1)-1,2:size(localI,1)-1);
        quiver( x, y ,  Ix(2:size(localI,1)-1,2:size(localI,1)-1), ...
            flipud(-Iy(2:size(localI,1)-1,2:size(localI,1)-1)) );
        axis([1, size(localI,1), 1, size(localI,1)])
        axis square
        title('gradients (with Gaussian weight)')
        subplot(2,2,4)
        plot(1:36,orientationHist, '*-');
        xlim([1, 36])
        ylim([0, 1.1*max(orientationHist)])
        %        plot(0:10:350,orientationHist);
        %        title(['orientation histogram,  \theta = ' num2str(keypoints(ct,4))/pi*180] );
        title('orientation histogram' );
        axis square
        waitforbuttonpress
    end
    
end

keypoints = keypoints(1:extraKeypointCt, :);
HessFailedKeypoints = HessFailedKeypoints(1:HessFailedCt, :);

%%
% Indicate the keypoints in the image.

%  I can make this a bit fancier by coloring the circles
%   'color',color_circle( idx )
Irgb = zeros( size(I, 1),  size(I, 2), 3);
Irgb(:,:,1) = I;
Irgb(:,:,2) = I;
Irgb(:,:,3) = I;
figure
imagesc( uint8(I) );
colormap gray
hold on;
title('Scale of keypoint is indicated by circle radius');
for ct = 1:extraKeypointCt
    viscircles([keypoints(ct,2) keypoints(ct,3)], ceil(sigmas(keypoints(ct,1))),'Color','r', 'LineWidth', 1, 'EnhanceVisibility', false);
    
    %  Draw the line indicating the dominant orientation
    line([keypoints(ct,2), keypoints(ct,2) + cos(keypoints(ct,4)) * sigmas(keypoints(ct,1))], ...
        [keypoints(ct,3), keypoints(ct,3) + sin(keypoints(ct,4)) *sigmas(keypoints(ct,1))],  'Color','r', 'LineWidth', 1)
end

for ct = 1:HessFailedCt
    viscircles([HessFailedKeypoints(ct,2) HessFailedKeypoints(ct,3)], ceil(sigmas(HessFailedKeypoints(ct,1))),'Color','b', 'LineWidth', 1, 'EnhanceVisibility', false);
    
    %line([keypoints(ct,2), keypoints(ct,2) + cos(keypoints(ct,4)) * sigmas(keypoints(ct,1))], ...
    %    [keypoints(ct,3), keypoints(ct,3) + sin(keypoints(ct,4)) *sigmas(keypoints(ct,1))],  'Color','r', 'LineWidth', 1)
end

if (1)
    figure
    imagesc( uint8(Irgb) );
    hold on;
    title('Scale of keypoint is indicated by circle radius');
end

%%   Rescale the image and rerun the code.

%   I = imresize(I, 0.8);


%%  assume a scalespace with 16 slices;   plot the slices
%
function subplotScalespace( scalespace, figtitle )
for j = 1:4
    for k = 1:4
        subplot(4,4, (j-1)*4+k)
        imagesc(squeeze(scalespace((j-1)*4+k,:,:)));
        colorbar
    end
end
sgtitle(figtitle)
%  trick from web   (not working quite right yet, ....  TODO)
%hp4 = get(subplot(4,4,16),'Position');
%colorbar('Position', [hp4(1)+hp4(3)+0.01  hp4(2)  0.1  hp4(2)+hp4(3)*2.1]);
end