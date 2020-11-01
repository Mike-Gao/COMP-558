% Load Image %
im_syn = rgb2gray(imread('synthetic.jpg'));
im_nat = rgb2gray(imread('mtl.jpg'));

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


% Q4: local extrema
figure
getpts(i)


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
lim = 50;
pts = zeros(100000,4);
cur = 1;
for i = 2:15
    sig = getsig(i);
    offset = ceil(2*sig)+1;
    for x = offset: w - offset
        for y = offset: h - offset
               if abs(dog(y,x,i)) > lim
                   neighbors = dog(y-1:y+1,x-1:x+1,i-1,i+1);
                   max = max(neighbors, [], 'all');
                   min = min(neighbors, [], 'all');
                   if (dog(y,x,i) == max || dog(y,x,i) == min)
                       pts(cur,1:3) = [x y sig];
                       cur = cur + 1;
                   end
               end
        end
    end
end

function plotimage(im,map,t)
    imagesc(im)
    colormap(map)
    colorbar
    title(t)
end
    
